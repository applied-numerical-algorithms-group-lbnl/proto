#include "Proto_Brick.H"
#include "Proto_Point.H"
#include <Python.h>
#include <cxxabi.h>
#include <iostream>
#include <memory>

void StencilProgram::compile() {
  std::stringstream command;

  command << BRICK_COMPILER << " " << BRICK_FLAGS << " " << name << ".cpp -o " << name
          << ".so -shared -fPIC" << std::flush;

  int ret = system(command.str().c_str());
  if (WEXITSTATUS(ret) != EXIT_SUCCESS) {
    std::cout << "compilation failed" << std::endl;
    exit(EXIT_FAILURE);
  }
}

bool StencilProgram::load() {
  dynlib = dlmopen(LM_ID_NEWLM, (name + ".so").c_str(), RTLD_NOW | RTLD_LOCAL);
  if (!dynlib)
    return false;

  const char *dlsym_error;

  fun_ptr = (compiled_t)dlsym(dynlib, "stencil");
  dlsym_error = dlerror();
  if (dlsym_error != nullptr)
    return false;

  auto desc_ptr = (char **)dlsym(dynlib, "stencil_desc");
  dlsym_error = dlerror();
  if (dlsym_error != nullptr)
    return false;
  stencil_desc = *desc_ptr;

  return true;
}

struct CompileRuntimeImplementation {
  PyObject *global;
  PyObject *local;
};

CompileRuntime::CompileRuntime() {
  Py_Initialize();

  implementation = new CompileRuntimeImplementation();

  auto typed_implem = (CompileRuntimeImplementation *)implementation;
  typed_implem->global = PyDict_New();
  typed_implem->local = PyDict_New();

  PyRun_String("import sys", Py_file_input, typed_implem->global, typed_implem->local);
  PyRun_String("print(sys.version)", Py_file_input, typed_implem->global, typed_implem->local);
}

void CompileRuntime::run(const std::string &code) {
  auto typed_implem = (CompileRuntimeImplementation *)implementation;
  PyRun_String(code.c_str(), Py_file_input, typed_implem->global, typed_implem->local);
}

std::string CompileRuntime::eval(const std::string &code) {
  auto typed_implem = (CompileRuntimeImplementation *)implementation;
  char *cstr;
  PyObject *pret =
      PyRun_String(code.c_str(), Py_eval_input, typed_implem->global, typed_implem->local);
  PyArg_Parse(pret, "s", &cstr);
  std::string ret = cstr;
  Py_DECREF(pret);
  return cstr;
}

CompileRuntime::~CompileRuntime() {
  Py_Finalize();
  delete (CompileRuntimeImplementation *)implementation;
};

void CompilationBase::loadCompiledFiles() {
  size_t pos = dir_prefix.rfind('/');

  // It must include directory "./" for current directory
  assert(pos != std::string::npos);

  std::string to_match = dir_prefix.substr(pos + 1);

  std::string dir_name = dir_prefix.substr(0, pos);

  DIR *dir = opendir(dir_name.c_str());
  if (dir != nullptr) {
    struct dirent *ent;
    while ((ent = readdir(dir)) != nullptr) {
      std::string fname(ent->d_name);
      if (fname == "." || fname == "..")
        continue;

      if (to_match == fname.substr(0, to_match.length())) {
        auto len = fname.length();
        if (len > 3 && fname.substr(len - 3) == ".so") {
          // This is a shared library
          std::string library_name = dir_prefix + fname.substr(0, len - 3);
          auto program = std::make_shared<StencilProgram>(library_name);
          size_t sz;
          try {
            int id = std::stoi(fname.substr(0, len - 3), &sz);
            if (sz == len - 3)
              files.insert(id);
            // Add to data base when loading is successful
            if (program->load()) {
              std::cout << "Loaded stencil " << program->getDesc() << std::endl;
              base[program->getDesc()] = program;
            }
          } catch (std::invalid_argument &e) {
          }
        }
      }
    }
    closedir(dir);
  }
}

void CompilationBase::initCompileRuntime() {
  compileRuntime = new CompileRuntime();

  // Python script to initialize some environment for latter use
  std::stringstream sstr;

  // The following will import the code generator into the python path
  sstr << "sys.path.append('" << BRICK_CODEGEN << "')" << std::endl;
  sstr << "from st.expr import Index, ConstRef, If, IntLiteral\nfrom st.grid import Grid\n"
       << std::endl;

  sstr << "from st.codegen.backend import Brick\n"
       << "from st.codegen.base import CodeGen\n"
       << "from io import StringIO\n"
       << std::endl;

  sstr << "def get_backend(vec):\n"
          "  from st.codegen.backend import BackendAVX512, BackendAVX2, BackendScalar, "
          "BackendCUDA, BackendFlex, BackendSSE, BackendCuFlex\n"
          "  if vec == 'AVX512':\n"
          "    return BackendAVX512()\n"
          "  elif vec == 'AVX2':\n"
          "    return BackendAVX2()\n"
          "  elif vec == 'SSE':\n"
          "    return BackendSSE()\n"
          "  elif vec == 'CUDA':\n"
          "    return BackendCUDA()\n"
          "  elif vec == 'FLEX':\n"
          "    return BackendFlex()\n"
          "  elif vec == 'CUFLEX':\n"
          "    return BackendCuFlex()\n"
          "  elif vec == 'OPENCL':\n"
          "    return BackendCUDA(16, 'sglid', ocl=True)\n"
          "  elif vec == 'SYCL':\n"
          "    return BackendCUDA(16, 'sglid', ocl=False)\n"
          "  elif vec == 'HIP':\n"
          "    return BackendCUDA(64, 'hipThreadIdx_x')\n"
          "  return BackendScalar()\n\n"
       << std::endl;
  compileRuntime->run(sstr.str());
}
