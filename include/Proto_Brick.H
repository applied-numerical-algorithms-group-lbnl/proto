#ifndef PROTO_PROTO_BRICK_H
#define PROTO_PROTO_BRICK_H

#ifdef PROTO_BRICK

#include "brick.h"
#include "bricksetup.h"
#include "immintrin.h"
#include "Proto_Point.H"
#include <dlfcn.h>
#include <unistd.h>
#include <utility>

#define BVEC "AVX2"

#if DIM == 3
#define BSIZE 8
#define BDIM BSIZE, BSIZE, BSIZE
#define BFOLD 2, 2
#elif DIM == 2
#define BSIZE 16
#define BDIM BSIZE, BSIZE
#define BFOLD 2, 2
#endif

typedef Brick<Dim<BDIM>, Dim<BFOLD>> Proto_brick;

namespace Proto {
template <unsigned dim, typename Func>
inline __attribute__((always_inline)) void
_iterate_partial_brick(const Point &a_p_start, const Point &a_p_end, const Func &a_point_f,
                       Point &a_p, StopTag t) {
  a_point_f(a_p);
}

template <unsigned dim, typename Func>
inline __attribute__((always_inline)) void
_iterate_partial_brick(const Point &a_p_start, const Point &a_p_end, const Func &a_point_f,
                       Point &a_p, RunningTag t) {
  auto d = dim - 1;
  for (a_p[d] = a_p_start[d]; a_p[d] < a_p_end[d]; ++a_p[d])
    _iterate_partial_brick<dim - 1>(a_p_start, a_p_end, a_point_f, a_p, TagSelect<dim - 1>::value);
}

template <unsigned dim, typename Func1, typename Func2>
inline __attribute__((always_inline)) void
_iterate_brick(const Point &a_lb, const Point &a_rb, const Func1 &a_point_f, const Func2 &a_brick_f,
               Point &a_p_start, Point &a_p_end, StopTag tag, StopTag t) {
  // Partial brick iteration
  Point p = a_p_start;
  _iterate_partial_brick<DIM>(a_p_start, a_p_end, a_point_f, p, TagSelect<DIM>::value);
}

template <unsigned dim, typename Func1, typename Func2>
inline __attribute__((always_inline)) void
_iterate_brick(const Point &a_lb, const Point &a_rb, const Func1 &a_point_f, const Func2 &a_brick_f,
               Point &a_p_start, Point &a_p_end, RunningTag tag, StopTag t) {
  a_brick_f(a_p_start);
}

template <unsigned dim, typename Func1, typename Func2, typename Tag>
inline // __attribute__((always_inline))
    void
    _iterate_brick(const Point &a_lb, const Point &a_rb, const Func1 &a_point_f,
                   const Func2 &a_brick_f, Point &a_p_start, Point &a_p_end, Tag tag,
                   RunningTag t) {
  // Rely on two's complement negative integers
  constexpr int mask_r = -BSIZE;

  auto d = dim - 1;
  // First peel
  a_p_start[d] = a_lb[d];
  a_p_end[d] = (a_lb[d] + BSIZE - 1) & mask_r;
  if (a_p_start[d] != a_p_end[d])
    _iterate_brick<dim - 1>(a_lb, a_rb, a_point_f, a_brick_f, a_p_start, a_p_end,
                            TagSelect<0>::value, TagSelect<dim - 1>::value);

  // Loop
  for (a_p_start[d] = a_p_end[d]; a_p_start[d] < (a_rb[d] & mask_r); a_p_start[d] += BSIZE) {
    a_p_end[d] = a_p_start[d] + BSIZE;
    _iterate_brick<dim - 1>(a_lb, a_rb, a_point_f, a_brick_f, a_p_start, a_p_end, tag,
                            TagSelect<dim - 1>::value);
  }
  // Last peel
  if (a_p_start[d] < a_rb[d]) {
    a_p_end[d] = a_rb[d];
    _iterate_brick<dim - 1>(a_lb, a_rb, a_point_f, a_brick_f, a_p_start, a_p_end,
                            TagSelect<0>::value, TagSelect<dim - 1>::value);
  }
}

/**
  \brief Apply function to all points in a range (rectangular) under bricks
  \tparam Func1 point-wise function type
  \tparam Func2 brick-wise function type
  \param lb lower bound (inclusive)
  \param rb upper bound (inclusive)
  \param point_f point-wise function
  \param brick_f brick-wise function
 */
template <typename Func1, typename Func2>
inline void iterate_brick(const Point &a_lb, const Point &a_rb, const Func1 &a_point_f,
                          const Func2 &a_brick_f) {
  Point p_start = Point::Zeros();
  Point p_end = Point::Zeros();
  _iterate_brick<DIM>(a_lb, a_rb, a_point_f, a_brick_f, p_start, p_end, TagSelect<1>::value,
                      TagSelect<DIM>::value);
}
} // namespace Proto

#ifdef PROTO_BRICK_JIT

#include <dirent.h>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <sys/stat.h>
#include <typeinfo>
#include <vector>

/// A stencil kernel program
class StencilProgram {
public:
  /// Type of the compiled program
  typedef void (*compiled_t)(Proto_brick &a_src, Proto_brick &a_dest, int brick_idx,
                             bool initToZero, std::vector<bElem> &coefs);

  /// Initialize a stencil kernel program using a string file name (without .so/.cpp)
  explicit StencilProgram(std::string name) : name(std::move(name)) {}

  /// Compile the kernel program (name + ".cpp")
  void compile();

  /// Load the kernel program with (name + ".so")
  bool load();

  /// Return the encoded stencil description. This is set from the loaded stencil program
  inline std::string getDesc() { return stencil_desc; }

  /**
   * @brief Execute this stencil program
   *
   * This is only a wrapper showing the different parameters to this problem. It is recommended to
   * use getFunction() to obtain the function pointer.
   *
   * @param a_src Source brick
   * @param a_dest Destination brick
   * @param brick_idx Index to the brick
   * @param initToZero When true, initialize the brick value to zero
   * @param coefs Coefficients to the stencil
   */
  inline void execute(Proto_brick &a_src, Proto_brick &a_dest, int brick_idx, bool initToZero,
                      std::vector<bElem> &coefs) {
    fun_ptr(a_src, a_dest, brick_idx, initToZero, coefs);
  }

  ~StencilProgram() {
    if (dynlib)
      dlclose(dynlib);
  }

  /**
   * @brief Get the Function pointer
   * @return void* Function pointer, type is ::compiled_t
   */
  inline void *getFunction() { return (void *)fun_ptr; }

private:
  void *dynlib = nullptr;       //!< the loaded shared library
  compiled_t fun_ptr = nullptr; //!< Function pointer
  std::string name;             //!< basename of stencil files
  std::string stencil_desc;     //!< encoded offset description of the stencil
};


/// Wrapper to the code generator
class CompileRuntime {
public:
  CompileRuntime();

  /**
   * Evaluating a piece of code and return the result in string
   * @param code
   * @return
   */
  std::string eval(const std::string &code);

  /**
   * Run a piece of code with no return value
   * @param code
   */
  void run(const std::string &code);

  ~CompileRuntime();

private:
  void *implementation;
};

/**
 * Interface to the compiled collection of stencils, new stencils will be automatically compiled
 * This object will control how each stencil is compiled and how they are stored
 */
class CompilationBase {
private:
  constexpr static const char *const encodeString = "0123456789abcdefghijklmnopqrstuvwxyz";
  constexpr static const int encodeOffset = 18;

  void loadCompiledFiles();

public:
  /**
   * Initialization
   * @param dir_prefix
   */
  explicit CompilationBase(std::string &&dir_prefix = "./")
      : dir_prefix(dir_prefix), compileRuntime(nullptr) {
    // Load all stencils
    loadCompiledFiles();
  }

  template <typename T> static std::string encodeName(const std::vector<Proto::Point> &points) {
    std::stringstream ret;
    ret << DIM << typeid(T).name();
    for (auto &pt : points)
      for (int i = 0; i < DIM; ++i)
        ret << encodeString[encodeOffset + pt.m_tuple[i]];
    return ret.str();
  }

  /**
   * Decoding the list of points from the name of the code
   * @param str the name of the code
   * @return a list of points
   */
  template <typename T> static std::vector<Proto::Point> decodeName(std::string &&str) {
    // Dimension must match
    assert(DIM == (str[0] - '0'));

    // The second character is the type that also should match
    assert(str[1] == typeid(T).name()[0]);

    // The points to return
    std::vector<Proto::Point> pts;

    for (int i = 2; i < str.length(); i += DIM) {
      int dims[DIM];
      for (int j = 0; j < DIM; ++j) {
        char c = str[i + j];
        if (c <= '9' && c >= '0')
          dims[j] = c - '0';
        else
          dims[j] = c - 'a' + 10;
        dims[j] -= encodeOffset;
      }
      pts.template emplace_back(dims);
    }

    return pts;
  }

  void initCompileRuntime();

  /**
   * @brief Get the StencilProgram
   * 
   * @tparam T stencil's program
   * @param points Stencil's offsets
   * @return std::shared_ptr<StencilProgram> 
   */
  template <typename T>
  std::shared_ptr<StencilProgram> getStencil(const std::vector<Proto::Point> &points) {
    auto basename = encodeName<T>(points);
    auto result = base.find(basename);
    if (result != base.end())
      return result->second;

    if (compileRuntime == nullptr) {
      initCompileRuntime();
    }

    int fileId = 0;
    if (files.size() != 0)
      fileId = std::max(fileId, *(files.rbegin()) + 1);
    std::string filename = dir_prefix + std::to_string(fileId);

    // Create the program and store in <name>.cpp
    // Python script to generate code
    std::stringstream sstr;

    // Declare indices
    for (int i = 0; i < DIM; ++i)
      sstr << (char)('i' + i) << " = Index(" << i << ")\n";

    sstr << std::endl;

    sstr << "input = Grid('src', " << DIM << ")\noutput = Grid('dest', " << DIM << ")\n"
         << std::endl;

    sstr << "coefs = []\nfor c in range(" << points.size()
         << "):\n  coefs.append(ConstRef('coefs[{}]'.format(c)))\n"
         << std::endl;

    sstr << "zero = ConstRef('initToZero')\n" << std::endl;

    for (int i = 0; i < points.size(); ++i) {
      if (i == 0)
        sstr << "calc = ";
      else
        sstr << " + ";
      sstr << "coefs[" << i << "] * input(";
      for (int j = 0; j < DIM; ++j) {
        if (j != 0)
          sstr << ", ";
        sstr << (char)('i' + j) << " + " << points[i].m_tuple[j];
      }
      sstr << ")";
    }

    sstr << std::endl;

    sstr << "output(";
    for (int i = 0; i < DIM; ++i) {
      if (i != 0)
        sstr << ", ";
      sstr << (char)('i' + i);
    }
    sstr << ").assign(If(zero, IntLiteral(0), output(";
    for (int i = 0; i < DIM; ++i) {
      if (i != 0)
        sstr << ", ";
      sstr << (char)('i' + i);
    }
    sstr << ")) + calc)\n" << std::endl;

    int prec = (typeid(T) == typeid(double)) ? 2 : 1;

    // Code generation options
    sstr << "vec = '" << BVEC << "'\n"
         << "backend = get_backend(vec)\n"
         << "dim = list(reversed(" << VS_TOSTR((BDIM,)) << "))\n"
         << "fold = list(reversed(" << VS_TOSTR((BFOLD,)) << "))\n"
         << "stride = [1] * len(dim)\n"
            "layout = Brick(dim=dim, fold=fold, prec="
         << prec << ", brick_idx='b',\n"
         << "               cstruct=False if vec != 'OPENCL' and vec != 'SYCL' else True)\n"
            "cg = CodeGen(backend, layout, dag_msize=5, scatter_thres=1.5, unroll=False,\n"
            "             dimsplit=True, stride=stride, K=2, klimit=20)\n"
            "s = StringIO()\n"
            "cg.gencode([output], s)\n"
         << std::endl;

    // Run the python code gen
    compileRuntime->run(sstr.str());

    std::string code = compileRuntime->eval("s.getvalue()");

    // Now put it into a "C++" program

    std::ofstream cppcode((filename + ".cpp").c_str());

    // The cppcode needs basic brick definitions and main function
    cppcode << "#define bElem double" << std::endl;
    cppcode << "#include \"" << BRICK_PATH << "/include/minimal_brick.h\"" << std::endl;
    cppcode << "#include <immintrin.h>" << std::endl;
    cppcode << "#include <vector>" << std::endl;
    cppcode << "const char * stencil_desc = \"" << basename << "\";" << std::endl;
    cppcode << "typedef Brick<Dim<" << VS_TOSTR(BDIM) << ">, Dim<" << VS_TOSTR(BFOLD)
            << ">> Proto_brick;" << std::endl;

    cppcode << "extern \"C\"\n"
               "void stencil(Proto_brick &src,\n"
               "             Proto_brick &dest,\n"
               "             int b,\n"
               "             bool initToZero,\n"
               "             std::vector<double> &coefs) {\n"
            << code << "}\n"
            << std::endl;

    cppcode.close();

    auto program = std::make_shared<StencilProgram>(filename);

    program->compile();
    program->load();

    base[basename] = program;

    return program;
  }

  ~CompilationBase() { delete compileRuntime; }

  static CompilationBase *getInstance() {
    static CompilationBase *instance = nullptr;
    if (!instance) {
      std::cout << "Initializing compilation base" << std::endl;
      auto status = mkdir("./stencils", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (status != 0 && errno != EEXIST)
        std::cout << "Error creating compilation folder" << status << std::endl;
      instance = new CompilationBase("./stencils/");
    }
    return instance;
  }

private:
  std::string dir_prefix;
  std::map<std::string, std::shared_ptr<StencilProgram>> base;
  std::set<int> files;
  CompileRuntime *compileRuntime;
};

#endif // PROTO_BRICK_JIT

#endif // PROTO_BRICK

#endif // PROTO_PROTO_BRICK_H
