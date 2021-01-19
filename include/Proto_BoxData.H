#ifndef _PROTO_RECTMDARRAY_H_
#define _PROTO_RECTMDARRAY_H_
#include "Proto_Timer.H"

#include "Proto_Box.H"
#include "Proto_PowerItoI.H"
#include "Proto_Stencil.H"
#include "Proto_Reduction.H"

#include <memory> //for shared_ptr
#include <iostream>
#include <iomanip> //for pretty printing
#include <limits> // for max / min functions
#include <cstdlib> // for size_t
#include <cmath> // floating point abs
#include <map>
#include <unordered_map>
#include "implem/Proto_Stack.H"
#include "Proto_DisjointBoxLayout.H"

// Please do not set this value to false (slow, race condition).
#define DEFAULT_USE_STACK true

namespace Proto
{
   enum BoxDataOp
  {
    Invalid = -1,
    Copy    = 0,
    Add     = 1,
    Assign,
    Subtract, 
    Multiply, 
    Divide, 
    NUM_OP_TYPES
  };
}

 
// This is to facilitate automated testing.
#ifdef PROTO_MEM_CHECK
class memcheck
{
public:
  static int numcopies;
  inline  static void FLUSH_CPY(){numcopies = 0;}
};
#endif
using std::shared_ptr;
namespace Proto 
{

// forward declarations
  template<typename T> class Stencil;
  template<typename T, unsigned int C, unsigned char D, unsigned char E>
  class LazyStencil;
  template <typename T, unsigned int C, unsigned char D, unsigned char E>
  class BoxData;

  inline void null_deleter_boxdata(void* ptr)
  {
  }

  inline void stack_deleter_boxdata(void* ptr)
   {
     StackAlloc& allocator  = DisjointBoxLayout::getStackAllocator();
     allocator.free(ptr);
   }


/** @defgroup pointwise_operations Pointwise Operations*/

/** @addtogroup pointwise_operations*/
/*@{*/

/// Component-Space Interval
/**
    Defines a subspace of tensor components.
    Used chiefly for copying and and linearized buffer IO.
    Usually the user will not need to explicitly create a CInterval.
    See the documentation for BoxData::copyTo for an example. 
*/
class CInterval 
{
public:
    /// Bounds Constructor
    /**
        Builds the interval:
          {{a_c0, a_c1},{a_d0, a_d1},{a_e0, a_e1}}
    */
    inline CInterval(unsigned int a_c0, unsigned int a_c1,
                     unsigned int a_d0 = 0, unsigned int a_d1 = 0,
                     unsigned int a_e0 = 0, unsigned int a_e1 = 0)
    {
      PROTO_ASSERT(a_c1 >= a_c0,
      "CInterval(...) invalid for a_c0 = %i and a_c1 = %i. \
      lower bound must be less than or equal to high bound.",
      a_c0, a_c1);
      PROTO_ASSERT(a_d1 >= a_d0,
      "CInterval(...) invalid for a_d0 = %i and a_d1 = %i. \
      lower bound must be less than or equal to high bound.",
      a_d0, a_d1);
      PROTO_ASSERT(a_e1 >= a_e0,
      "CInterval(...) invalid for a_e0 = %i and a_e1 = %i. \
      lower bound must be less than or equal to high bound.",
      a_e0, a_e1);
      m_comps[0][0] = a_c0; m_comps[0][1] = a_c1;
      m_comps[1][0] = a_d0; m_comps[1][1] = a_d1;
      m_comps[2][0] = a_e0; m_comps[2][1] = a_e1;
    }
    /// List Constructor
    /**
        Build a CInterval using the syntax:
        @code
        CInterval I0{{c0, c1},{d0, d1},{e0, e1}};
        // OR
        CInterval I1({{c0, c1},{d0, d1},{e0, e1}});
        // Or, for a single component index:
        CInterval I2{c0, c1};
        // OR
        CInterval I3({c0, c1});
        @endcode
    */
    CInterval(std::initializer_list<
             std::initializer_list<unsigned int>> a_lst)
    {
      int comp = 0;
      for (auto& c : a_lst)
      {
        if (c.size() == 0)
        {
          m_comps[comp][0] = 0;
          m_comps[comp][1] = 0;
          comp++;
          continue;
        }
        int bound = 0;
        for (auto& b : c)
        {
          m_comps[comp][bound] = b;
          bound++;
        }
        comp++;
      }
      PROTO_ASSERT(m_comps[0][1] >= m_comps[0][0],
      "CInterval(...) invalid for a_c0 = %i and a_c1 = %i. \
      lower bound must be less than or equal to high bound.",
      m_comps[0][0], m_comps[0][1]);
      PROTO_ASSERT(m_comps[1][1] >= m_comps[1][0],
      "CInterval(...) invalid for a_d0 = %i and a_d1 = %i. \
      lower bound must be less than or equal to high bound.",
      m_comps[1][0], m_comps[1][1]);
      PROTO_ASSERT(m_comps[2][1] >= m_comps[2][0],
      "CInterval(...) invalid for a_e0 = %i and a_e1 = %i. \
      lower bound must be less than or equal to high bound.",
      m_comps[2][0], m_comps[2][1]);
    }
    
    /// Lower Bound
    /**
        Retrieve the lower bound along the given component axis
    
        \param a_comp   A component axis in [0,3)
    */
    inline unsigned int low(unsigned int a_comp) const
    {
      PROTO_ASSERT((a_comp < 3),
      "CInterval::low(unsigned int a_comp) invalid for a_comp = %i.\
      a_comp must be in [0,3).",
      a_comp);
      return m_comps[a_comp][0];
    };
    
    /// Upper Bound
    /**
        Retrieve the upper bound component
        for component axis <code> a_comp </code>
        
        \param a_comp   A component axis in [0,3)
    */
    inline unsigned int high(unsigned int a_comp) const
    {
      PROTO_ASSERT((a_comp < 3),
      "CInterval::high(unsigned int a_comp) invalid for a_comp = %i.\
      a_comp must be in [0,3).",
      a_comp);
      return m_comps[a_comp][1];
    };
    
    /// Contains Query
    /**
        Query if <code>*this</code> contains <code>index</code>
        of component <code>comp</code>
        
        \param a_index  An index
        \param a_comp   A component axis in [0,3)
    */
    inline bool contains(unsigned int a_index,
                         unsigned int a_comp) const
    {
        PROTO_ASSERT((a_comp < 3),
        "CInterval::contains(unsigned int a_comp) \
        invalid for a_comp = %i. a_comp must be in [0,3).",
        a_comp);
        return ((a_index >= m_comps[a_comp][0])
             && (a_index <= m_comps[a_comp][1]));
    }

    /// Size Query
    /**
        Returns the number of components in <code> *this </code>
        along the component axis <code> a_comp </code>
        
        \param a_comp   A component axis in [0,3)
    */
    inline unsigned int size(unsigned int a_comp) const
    {
        PROTO_ASSERT((a_comp < 3),
        "CInterval::size(unsigned int a_comp) \
        invalid for a_comp = %i. a_comp must be in [0,3).",
        a_comp);
        return (m_comps[a_comp][1] - m_comps[a_comp][0] + 1);
    }

    /// Print
    inline void print() const
    {
      std::cout << "{ ";
      std::cout << "{" << m_comps[0][0] << ", " << m_comps[0][1] << "}";
      std::cout << ", {" << m_comps[1][0] << ", " << m_comps[1][1] << "}";
      std::cout << ", {" << m_comps[2][0] << ", " << m_comps[2][1] << "}";
      std::cout << " }";
    }
private:
    unsigned int m_comps[3][2]; ///< bounds of the interval
};

/// CInterval IOStream Operator
inline std::ostream& operator<<(std::ostream& a_os,
                                const CInterval& a_int)
{
      a_os << "{ ";
      a_os << "{" << a_int.low(0) << ", " << a_int.high(0) << "}";
      a_os << ", {" << a_int.low(1) << ", " << a_int.high(1) << "}";
      a_os << ", {" << a_int.low(2) << ", " << a_int.high(2) << "}";
      a_os << " }";
      return a_os;
}

/// Pointwise Variable
/**
    Var is used for implementing pointwise operators on BoxDatas.
    See the documentation for forall for more information.
    Template parameters refer to the type and dimension of data held by the associated BoxData.
    See BoxData documentation.
    
    \tparam T Type of data in array (bool, int, double, etc.)
    \tparam C Number of components in first data index. Defaults to 1. 
    \tparam D Number of components in second data index. Defaults to 1. 
    \tparam E Number of components in third data index. Defaults to 1.
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1>
  class Var
  {
  public:

    /// Pointwise Accessor
    /**
       Access component (c,d,e) of the <code>BoxData<T,C,D,E></code> associated with <code>*this</code>.
       
       \param a_c   First component index
       \param a_d   Second component index (default: 0)
       \param a_e   Third component index  (default: 0)
    */
    CUDA_DECORATION
    inline T& operator()(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0)  
    {
#ifdef __CUDA_ARCH__
      int idx = threadIdx.x + blockIdx.x*blockDim.x;
      int idy = threadIdx.y + blockIdx.y*blockDim.y;
#if DIM == 3
      int idz = threadIdx.z + blockIdx.z*blockDim.z;
      return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + boxDimX * (idy + idz * boxDimY)];
#else
      return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + idy * boxDimX];
#endif
#else
      return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
#endif
    }

    /// Pointwise Accessor (Const)
    /**
       Access component (c,d,e) of the <code>const BoxData<T,C,D,E></code> associated with <code>*this</code>.
       
       \param a_c   First component index
       \param a_d   Second component index (default: 0)
       \param a_e   Third component index  (default: 0)
    */
    CUDA_DECORATION
    inline const T& operator()(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0) const 
    {
#ifdef __CUDA_ARCH__
      int idx = threadIdx.x + blockIdx.x*blockDim.x;
      int idy = threadIdx.y + blockIdx.y*blockDim.y;
#if DIM == 3
      int idz = threadIdx.z + blockIdx.z*blockDim.z;
      return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + boxDimX * (idy + idz * boxDimY)];
#else
      return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx + idy * boxDimX];
#endif
#else
      return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
#endif
    }
 
    CUDA_DECORATION
    inline Var&  operator+=(unsigned int& a_increment) 
    {
      for (int ii = 0; ii < C*D*E; ii++)
      {
        m_ptrs[ii] += a_increment;
      }
      return *this;
    }

    // test
    CUDA_DECORATION
    inline Var&  operator+=(Point& a_p)
    {
#if DIM == 3
      unsigned int shift = a_p[0] + (a_p[1] + a_p[2]*boxDimY)*boxDimX;
#else
      unsigned int shift = a_p[0] + a_p[1]*boxDimX;
#endif
      return *this+=shift;
    }

    CUDA_DECORATION
    inline Var&  operator+=(const Point& a_p)
    {
#if DIM == 3
      unsigned int shift = a_p[0] + (a_p[1] + a_p[2]*boxDimY)*boxDimX;
#else
      unsigned int shift = a_p[0] + a_p[1]*boxDimX;
#endif
      return *this+=shift;
    }
   
    CUDA_DECORATION
    inline Var& operator++() 
    {
      for (int ii = 0; ii < C*D*E; ii++)
      {
        ++m_ptrs[ii];
      }
      return *this;
    }

    CUDA_DECORATION
    inline Var& operator--() 
    {
      for (int ii = 0; ii < C*D*E; ii++)
      {
        --m_ptrs[ii];
      }
      return *this;
    }

    unsigned int boxDimX;
    unsigned int boxDimY;
    unsigned int subBoxDimX;
    unsigned int subBoxDimY;
    T* m_ptrs[C*D*E];
  };

/// Multidimensional Rectangular Array
/**
    BoxData is the main dataholder class of Proto. It contains an 
    array of data of type T on a domain defined by a Box.
    Each data point is associated with a Point in the Box.
    C, D, and E define the tensor structure of the data:
    C,D,E = 1:            Data is scalar (default)
    C = N, D,E = 1:       Data is vector valued with N components
    C = M, D = N, E = 1:  Data is MxN matrix valued

    \tparam T Type of data in array (bool, int, double, etc.)
    \tparam C (Optional) Size of first component axis. Defaults to 1. 
    \tparam D (Optional) Size of second component axis. Defaults to 1. 
    \tparam E (Optional) Size of third component axis. Defaults to 1.
*/
  template <class T=double, unsigned int C=1, unsigned char D=1, unsigned char E=1>
  class BoxData 
  {
  public:

    typedef Var<T,C,D,E> reference;
    typedef BoxData<double,1,1,1> Scalar;
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Constructors and Define*/
    ///@{
  
    /// Copy constructor/assignment forbidden for all the usual reasons.
    BoxData(const BoxData<T,C,D,E>& a_src) = delete;
    BoxData& operator=(const BoxData<T,C,D,E>& a_src) = delete;

    /// Default Constructor
    /**
       Allocates no space; user must call define on resulting object.
    */
    BoxData();

    /// Box Constructor
    /**
      Creates and allocates an uninitialized array based on a Box input. 
      
      \param a_box              Box defining the domain of *this
      \param a_stackAllocation  If true, memory will be allocated with <code>Proto::Stack</code> (default: DEFAULT_USE_STACK); 
    */
    BoxData(const Box& a_box, bool a_stackAllocation=DEFAULT_USE_STACK);
    
    /// Box Constructor With Initialization
    /**
       Identical to calling the Box constructor followed by setVal(a_init).
       Not recommended in performance-critical code if initialization is unnecessary.
        
       \param a_box    Box defining the domain of *this
       \param a_init   T value used to initialize all data values
    */
    BoxData(const Box& a_box, T a_init);
    

    /// Move constructor
    /**
      Necessary in cases where we return BoxData by value, but don't want an  actual deep copy
      
      \param a_src    Source data
    */
    BoxData(BoxData<T,C,D,E>&& a_src);
    
    /// Box Define
    /**
       Allocates or reallocates memory for a previously declared BoxData.
       Any existing data in <code> *this </code> is thrown away.
       
       \param a_box              Box defining the domain of *this
       \param a_stackAllocation  If true, memory will be allocated with <code>Proto::Stack</code> (default: DEFAULT_USE_STACK); 
    */  
    void define(const Box& a_box, bool a_stackAllocation=DEFAULT_USE_STACK);


    ///define using malloc instead of stack allocation
    /** @private */
    void defineMalloc(const Box& a_box);
   

    ///define using stack allocation
    /** @private */
    void defineStack(const Box& a_box);
    
    ///define for aliasing
    /**@private */
    template<unsigned int ncomp>
    void define(BoxData<T,  ncomp>& a_input,
                unsigned int      & a_comp)
    {
      T* dataptr = a_input.dataPtr(a_comp);
      m_data   = a_input.getData();
      m_box    = a_input.box();
      m_rawPtr = dataptr;
    }           
   

    /// Raw Pointer Constructor
    /**
       Builds a vector valued BoxData by aliasing a raw pointer to T which should be of size
       <code> a_box.size()*a_ncomp </code>. The data in <code> a_ptr </code> is not copied.
       This constructor is not recommended for public use, but is provided to facilitate
       compatability with third party libraries.
       
       \param a_ptr   Source buffer of type T to which *this will be aliased
       \param a_box   Box defining the domain of *this
       \param a_ncomp (Optional) Dummy input for the number of components. Completely unused.
    */
    BoxData(const T* a_ptr, const Box& a_box, int a_ncomp = C);


    /// Raw Pointer Define
    /**
       Defines a vector valued BoxData by aliasing a raw pointer to T which should be of size
       <code> a_box.size()*a_ncomp </code>. The data in <code> a_ptr </code> is not copied.
       This constructor is not recommended for public use, but is provided to facilitate
       compatability with third party libraries.
       
       \param a_ptr   Source buffer of type T to which *this will be aliased
       \param a_box   Box defining the domain of *this
       \param a_ncomp (Optional) Dummy input for the number of components. Completely unused.
    */
    void define(const T* a_ptr, const Box& a_box, int a_ncomp = C);


    /// LazyStencil Constructor
    /**
        @private
        Builds a BoxData from a lazily evaluated Stencil operation. Called implicitly.
    
        \param a_op   A LazyStencil created by Proto::Stencil::operator() or Proto::InterpStencil::operator()
    */
    BoxData(LazyStencil<T,C,D,E>&& a_op);
   

    /// Slice Constructor
    /**
        @private
        Used internally. Not part of the public interface.
    */
    BoxData(shared_ptr<T> a_data,const T* a_ptr,const Box& a_box);


    /// Destructor.
    ~BoxData(); 
  
    ///@}
   
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Data Movement */
    ///@{

    template<BoxDataOp>
    void operatorT(const BoxData<T,C,D,E>& a_src);

    template<BoxDataOp>
    void operatorT(const T a_scale);

    /// Move Assignment Operator
    /**
       Moves data in a_src to *this without performing a deep copy.
    
       \param a_src   Source Data 
    */
    BoxData& operator=(BoxData<T,C,D,E>&& a_src);

    /// Copy on Intersection.
    /**
       Copy data from <code> *this </code> to <code> a_dest </code> 
       within the domain intersection. Does nothing if the domain intersection is empty.
       
       \param a_dest  Destination data holder
    */
    void copyTo(BoxData<T,C,D,E>& a_dest) const;

    /// Copy Region
    /**
       Copy with a prescribed Box argument and optional shift. Explicitly, this function copies 
       the subset of data from <code> *this </code> contained in <code> a_srcBox </code>
       into the region of <code> a_dest </code> defined by <code> a_srcBox.shift(a_destShift) </code>.

       This function will fail if <code> a_srcBox </code> is not contained in both <code> this->box() </code>
       and <code> a_dest.box().shift(a_destShift) </code>.

       \param a_dest       Destination data holder
       \param a_srcBox     Region of data to copy from *this
       \param a_destShift  (Optional) Determines region of a_dest to copy data to.
    */
    void copyTo(BoxData<T,C,D,E>&   a_dest,
                const Box&           a_srcBox,
                const Point&        a_destShift = Point::Zeros()) const;
    
    ///Copy with src and dest boxes.                                                                                                                                                                                            
    /**
        Copies from a <code> a_srcBox </code> in <code> *this </code> into 
        a <code> a_destBox </code> in <code> a_dest </code>. Used in ::Proto::Copier.                                                                                          
    */
    void copyTo(BoxData<T,C,D,E>& a_dest,
                const Box&    a_srcBox,
                const Box&    a_destBox) const;

    /// General Copy
    /**
        The most general form of copyTo.
        Copies a prescribed set of components from <code> *this </code>
        into <code> a_dest </code> within prescribed region with a possible shift. 
       
        Example Usage:
        \snippet Snippets.cpp proto_copyTo 
  
        \param a_dest       Destination data holder
        \param a_srcBox     Region of data to copy from *this
        \param a_srcComps   Components of *this to copy from                                                            
        \param a_destShift  Determines region of a_dest to copy data to (e.g. <code>a_srcBox.shift(a_destShift)</code>)
        \param a_destComps  Components of a_dest to copy into. Must be the same size as a_srcComps.
    */
    template <unsigned int Cdest,
              unsigned char Ddest,
              unsigned char Edest>
    void copyTo(BoxData<T,Cdest,Ddest,Edest>& a_dest,
                const Box&    a_srcBox,
                CInterval    a_srcComps,
                const Point& a_destShift,
                CInterval    a_destComps) const;

    /// General Copy From
    /**
        Nearly identical to the most general version of <code> copyTo </code> except the copy direction is reversed.
        Only valid for vector valued BoxData. This function is provided to facilitate interoperability with third party libraries, and is not recommended for public use unless necessary.
    
        \param a_dsrc     Source data
        \param a_srcBox   Source domain to copy from
        \param a_srcComp  First component to copy from source
        \param a_destBox  Destination domain to copy into
        \param a_destComp First component to copy into destiation
        \param a_numcomp  Number of components to copy
    */
    template<unsigned int Csrc>
    void copy(const BoxData<T,Csrc,D,E>& a_dsrc,
              const Box&                  a_srcBox,
              unsigned int               a_srcComp,
              const Box&                  a_destBox,
              unsigned int               a_destComp,
              unsigned int               a_numcomp);
    
    ///@}
  
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Accessors */
    ///@{
    
    /// Read Only Data Accessor 
    /**
      Read-only access to data stored in <code> *this </code>.
      This function is read-only to facilitate cross-platform code (e.g. GPU compatability).
      This function is for debugging ONLY, and will not work on the device.
      Pointwise write operations should be done through forall.

      \param a_pt   A Point in the domain of <code> *this </code>
      \param a_c    Value of first tensor index
      \param a_d    Value of second tensor index
      \param a_e    Value of third tensor index
    */
    CUDA_DECORATION
    inline T& operator()(const Point& a_pt,
                        unsigned int  a_c = 0,
                        unsigned char a_d = 0,
                        unsigned char a_e = 0) const
    {
      PROTO_ASSERT(m_box.contains(a_pt),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_pt not in this->box()");
      PROTO_ASSERT((a_c < C),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_c = %i. a_c must be in [0,%i)",a_c,C);
      PROTO_ASSERT((a_d < D),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_d = %i. a_d must be in [0,%i)",a_d,D);
      PROTO_ASSERT( (a_e < E),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_e = %i. a_e must be in [0,%i)",a_e,E);
      return m_rawPtr[index(a_pt,a_c,a_d,a_e)];
    }
  
    /// Index Accessor (Non-Const)
    /**
       Access *this using a linear index. Useful for interacting with
       BoxData like a regular buffer, but not recommended.

       \param a_index   Index in [0,this->size() = m_box.size()*C*D*E)
    */
    inline T* operator[](unsigned int a_index); 
    

    /// Index Accessor (Const)
    /**
        Access *this using a linear index. Useful for interacting with
        BoxData like a regular buffer, but not recommended.

        \param a_index   Index in [0,this->size() = m_box.size()*C*D*E)
    */
    inline const T* operator[](unsigned int a_index) const; 
   

    /// Read-Write Accessor (Const)
    /**
        Return a const pointer to a data point in <code> *this </code> given a Point and tensor indices.
        
        \param a_p  Point in the domain of this 
        \param a_c  (Optional) First tensor index. Defaults to 0.
        \param a_d  (Optional) Second tensor index. Defaults to 0.
        \param a_e  (Optional) Third tensor index. Defaults to 0.
    */
    inline const T* data(const Point& a_p,
                         unsigned int a_c = 0,
                         unsigned int a_d = 0,
                         unsigned int a_e = 0) const
    {
      size_t ioff = index(a_p, a_c, a_d, a_e);
      return &(m_rawPtr[ioff]);
    }

    /// Buffer Accessor (Non-Const)
    /** 
        Return the contiguous data buffer where the data in <code>*this</code> is stored.
        Not recommended for public use. 
    */
    inline T* data() 
    {
       return m_rawPtr;
    }

    /// Buffer Accessor (Const)
    /** 
        Return the contiguous data buffer where the data in <code>*this</code> is stored.
        Not recommended for public use. 
    */
    inline const T* data() const
    {
       return m_rawPtr;
    }

    /// Read-Write Accessor (Non-Const)
    /**
        Return a pointer to a data point in <code> *this </code> given a Point and tensor indices.
        
        \param a_p  Point in the domain of this 
        \param a_c  (Optional) First tensor index. Defaults to 0.
        \param a_d  (Optional) Second tensor index. Defaults to 0.
        \param a_e  (Optional) Third tensor index. Defaults to 0.
    */
    inline T* data(const Point& a_p,
                   unsigned int a_c = 0,
                   unsigned int a_d = 0,
                   unsigned int a_e = 0)
    {
      size_t ioff  = index(a_p, a_c, a_d, a_e);
      return &(m_rawPtr[ioff]);
    }
    
    /// Read-Write Accessor (Const)
    /**
        Return a pointer to the region of the data buffer where a given component starts.
        
        \param a_c  (Optional) First tensor index. Defaults to 0.
        \param a_d  (Optional) Second tensor index. Defaults to 0.
        \param a_e  (Optional) Third tensor index. Defaults to 0.
    */
    inline const T* dataPtr(unsigned int a_c = 0,
                            unsigned int a_d = 0,
                            unsigned int a_e = 0) const
    {
      Point pt = m_box.low();
      size_t ioff = index(pt, a_c, a_d, a_e);
      return &(m_rawPtr[ioff]);
    }

    /// Read-Write Accessor (Non-Const)
    /**
        Return a pointer to the region of the data buffer where a given component starts.
        
        \param a_c  (Optional) First tensor index. Defaults to 0.
        \param a_d  (Optional) Second tensor index. Defaults to 0.
        \param a_e  (Optional) Third tensor index. Defaults to 0.
    */
    inline T* dataPtr(unsigned int a_c = 0,
                      unsigned int a_d = 0,
                      unsigned int a_e = 0)
    {
      Point pt = m_box.low();
      size_t ioff = index(pt, a_c, a_d, a_e);
      return &(m_rawPtr[ioff]);
    }

    /// Access Shared Pointer
    /**
        Used internally. Not recommended for public use. 
    */  
    inline shared_ptr<T> getData() const
    {
      return m_data; ///< Data array
    }

    /// Compute Index
    /**
        Computes the index in <code>[0,this->size())</code> associated with a given Point and components
      
        \param a_pt   A Point in the domain of <code> *this </code>
        \param a_c    Value of first tensor index
        \param a_d    Value of second tensor index
        \param a_e    Value of third tensor index
    */
    CUDA_DECORATION
    inline size_t index(const Point a_pt,
                        unsigned int a_c = 0,
                        unsigned int a_d = 0,
                        unsigned int a_e = 0) const
    {
      PROTO_ASSERT(m_box.contains(a_pt),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
                   invalid for a_pt not in this->box()");
      PROTO_ASSERT((a_c < C),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
                    invalid for a_c = %i. a_c must be in [0,%i)",a_c,C);
      PROTO_ASSERT((a_d < D),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
                    invalid for a_d = %i. a_d must be in [0,%i)",a_d,D);
      PROTO_ASSERT((a_e < E),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
                    invalid for a_e = %i. a_e must be in [0,%i)",a_e,E);
      size_t m = m_box.size();
      size_t k = m_box.index(a_pt);
      return k+m*a_c+a_d*m*C+a_e*m*C*D;
    }
  

    // Get Domain Box
    /**
       Return's the domain Box of <code> *this </code>. 
       Includes ghost cells if it was built it that way.
    */
    inline Box box() const {return m_box;};
    
    // Get Size
    /**
        Returns the number of data points in <code> *this </code>.
        Return value is equal to <code> m_box.size()*C*D*E </code>
    */

    inline std::size_t size() const {return m_box.size()*C*D*E;};

    /// Defined Query
    /**
        Returns true if this is defined (e.g., it has memory allocated to it).
        This will return false for default constructed objects and true otherwise
    */
    inline bool defined() const {return bool(m_data);};
   

    /// Create Pointwise Access Variable (Non-Const)
    /**
       Not for public use
    */
    inline Var<T,C,D,E> var(const Point& a_pt);
    
    /// Create Pointwise Access Variable (Const)
    /**
       Not for public use
    */
    inline Var<T,C,D,E> var(const Point& a_pt) const;
 
    ///@}
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Algebraic Operations */
    ///@{

    /// Pointwise Addition on Intersection
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_rhs    Another BoxData
    */
    BoxData<T,C,D,E>& operator+=(const BoxData<T,C,D,E>& a_rhs);
    
    /// Pointwise Subtraction on Intersection
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_rhs    Another BoxData
    */
    BoxData<T,C,D,E>& operator-=(const BoxData<T,C,D,E>& a_rhs);

    /// Pointwise Multiplication on Intersection
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_rhs    Another BoxData
    */
    BoxData<T,C,D,E>& operator*=(const BoxData<T,C,D,E>& a_rhs);

    /// Pointwise Division on Intersection
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_rhs    Another BoxData
    */
    BoxData<T,C,D,E>& operator/=(const BoxData<T,C,D,E>& a_rhs);

    /// Pointwise Addition by Scalar
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_scale  A scalare of type T
    */
    BoxData<T,C,D,E>& operator+=(T a_scale);

    /// Pointwise Subtraction by Scalar
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_scale  A scalare of type T
    */
    BoxData<T,C,D,E>& operator-=(T scale);

    /// Pointwise Multiplication by Scalar
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_scale  A scalare of type T
    */
    BoxData<T,C,D,E>& operator*=(T scale);

    /// Pointwise Division by Scalar
    /**
        Not recommended for performance unless necessary.
        Try to fuse arithmetic operations into forall or Stencil operations if possible.
    
        \param a_scale  A scalare of type T
    */
    BoxData<T,C,D,E>& operator/=(T scale);
  
    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Utility */
    ///@{

    /// Initialize All Values
    /**
       Avoid calling this function as much as possible in performance critical code

       \param a_val A constant value. 
    */
    void setVal(const T& a_val);
  
    /// Set All Values in Box
    /**
       Avoid calling this function as much as possible in performance critical code
        
       \param a_val    A constant value.
       \param a_box    Domain to set to a_val
    */
    void setVal(const T&    a_val,
                const Box& a_box);

    /// Set All Values of Component in Box
    /**
       Avoid calling this function as much as possible in performance critical code
        
       \param a_val    A constant value. 
       \param a_box    Domain to set to a_val
       \param a_c      First index of component to set to a_val.
       \param a_d      Second index of component to set to a_val.
       \param a_e      Third index of component to set to a_val.
    */
    void setVal(const T& a_val,
                const Box& a_box,
                int a_c,
                int a_d = 0,
                int a_e = 0);

    /// Absolute Maximum Value (Global)
    /**
    GPU: calculates the maximum absolute value over the entire data set
    CPU: calls absMax() and updates a_Rxn.host if larger
    */
    void absMax(Reduction<T,Abs>& a_Rxn) const;

    /// Absolute Maximum Value (Componentwise)
    /**
       Returns the maximum absolute value of a given component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T absMax(int a_c = 0, int a_d = 0, int a_e = 0) const;

    /// Minimum Value (Componentwise)
    /**
       Returns the minimum value of a given component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T min(int a_c = 0, int a_d = 0, int a_e = 0) const;

    /// Maximum Value (Componentwise)
    /**
       Returns the maximum value of a given component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T max(int a_c = 0, int a_d = 0, int a_e = 0) const;

    /// sum (Componentwise)
    /**
       Return the sum of a particular component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T sum(int a_c = 0, int a_d = 0, int a_e = 0) const;

    /// Shift Domain
    /**
       Shifts the domain Box of this, moving all data along with it.
       e.g. data associated with Point p will now be associated with Point p + a_pt.

       \param a_pt  A Point interpreted as a shift vector.
    */
    inline void shift(const Point a_pt)
    {
      m_box = m_box.shift(a_pt);
    };
  
    /// Buffer Write
    /**
       Write a subset of data from <code> *this </code> into a C-array buffer.
       
       Example Usage:
       \snippet Snippets.cpp proto_linearInOut
       
       \param a_buf    Destination buffer
       \param a_box    Domain to write from
       \param a_comps  Components to write from
    */
    void linearOut(void*     a_buf,
                   const Box& a_box,
                   CInterval a_comps) const;

    void linearOut(void*                  a_buf,
                   const ::Proto::Box&      a_bx,
                   unsigned int a_startcomp,
                   unsigned int a_numcomps) const;


    /// Buffer Read
    /**
       Read data into *this from a C-array buffer populated by BoxData<T,C,D,E>::linearOut(...).
       The Box and components may be shifted with respect to those used for the read operation.
       (The number of component indices along each axis and the shape and size of the Box must be the same).
       See the example below.

       Example Usage:
       \snippet proto_linearInOut

       \param a_buf    Source buffer
       \param a_box    Domain to read from
       \param a_comps  Components to read from
    */
    void linearIn(void*     a_buf,
                  const Box& a_box,
                  CInterval a_comps);

    void linearIn(void*                  a_buf,
                  const ::Proto::Box&      a_bx,
                  unsigned int a_startcomp,
                  unsigned int a_numcomps);
 

   
    /// Contains CInterval
    /**
        Checks if an CInterval object is a subset of the component space
        of *this. This function is mostly used internally, but it can
        be useful for checking inputs to CopyTo and the slicing functions
    */
    inline bool contains(CInterval a_interval) const
    {
      bool ret = true;
      ret &= (a_interval.high(0) < C);
      ret &= (a_interval.high(1) < D);
      ret &= (a_interval.high(2) < E);
      return ret;
    }
    
    /// Check Aliasing
    /**
        Returns true if <code> *this </code> and <code> a_src </code> are aliased
        to the same data buffer. This function also returns true if one array
        is a slice of another.
    */
    template<unsigned int CC, unsigned char DD, unsigned char EE>
    inline bool isAlias(const BoxData<T, CC, DD, EE>& a_src) const
    {
      return (m_data.get() == a_src.m_data.get());
    }

    /// Print
    /**
        Default Print method. Outputs Domain Box and extrema of *this
    */
    void print() const;

    /// Print Data
    /**
        Pretty prints *all* of the data in *this. Prettiness may vary depending
        on the domain size. Generally looks best with DIM = 2 and Box edge lengths less than ~16.
        This function is purely for debugging purposes and should not be used in performance code.

        \param a_prec   (Optional) Desired precision for fixed decimal output. Defaults to 2.
    */    
    void printData(int a_prec = 2) const;
    
    /// Print Data in Box
    /**
        Pretty prints the data in *this within a given Box for all components.
        Input Boxes of size 16^DIM or smaller are recommended to maintain prettiness.
        This function is purely for debugging purposes and should not be used in performance code.
        
        \param a_box    Desired subset for printing
        \param a_prec   (Optional) Desired precision for fixed decimal output. Defaults to 2.
    */    
    void printData(const Box& a_box, int a_prec = 2) const;

    /// Print Component Data in Box
    /**
        Pretty prints the data in *this within a given Box for a single component.
        Input Boxes of size 16^DIM or smaller are recommended to maintain prettiness.
        This function is purely for debugging purposes and should not be used in performance code.
    
      
        \param a_box    Desired subset for printing
        \param a_c      First tensor index. Defaults to 0.
        \param a_d      Second tensor index. Defaults to 0.
        \param a_e      Third tensor index. Defaults to 0.
        \param a_prec   (Optional) Desired precision for fixed decimal output. Defaults to 2
        
    */    
    void printData(const Box& a_box, int a_c, int a_d, int a_e, int a_prec = 2) const;

    ///@}

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Friends
    
    template<class TT>
    friend class Stencil;
    
    template<class TT, unsigned int CC, unsigned char DD, unsigned char EE>
    friend BoxData<TT,CC,DD,EE> alias(BoxData<TT,CC,DD,EE>& a_original,
                                      const Point&              a_shift);

    template<class TT, unsigned int CC, unsigned char DD, unsigned char EE>
    friend const BoxData<TT,CC,DD,EE> alias(const BoxData<TT,CC,DD,EE>& a_original,
                                            const Point&                a_shift);

   class Interval;
   class IntVect;
   size_t size(const Box&      a_box,
                      const size_t a_comps) const
   {
     return a_box.size() * (sizeof(T) * a_comps);
   }



   size_t charsize(const ::Proto::Box&      a_bx,
                          unsigned int a_startcomp,
                          unsigned int a_numcomps) const {

    Box chbox(a_bx.low(), a_bx.high());

    const size_t interv =  a_startcomp + a_numcomps + a_startcomp + 1;
    return size(chbox, interv);
  }

  static int preAllocatable()
  {
    return 0; // static preAllocatable
  }


  private:

    bool                     m_stackAlloc;
    Box                      m_box; ///< Box defining the domain of *this
    ::std::shared_ptr<T>    m_data; ///< Data array
    T*                      m_rawPtr; ///< Raw pointer to the data
  };

///////////////////////////////////////////////////////////////////////////////////////////////
/** @name Alias and Slice Operators 
    The alias and slice operations facilitate BoxData operations while avoiding unnecessary copies.
    See the sample code below for an explanation of the syntax.

    Example <code> alias </code> usage:
    \snippet Snippets.cpp proto_alias
    
    Example <code> slice </code> usage:
    \snippet Snippets.cpp proto_slice
*/

///@{

/// Alias (Non-Const)
/**
   Creates a read-write alias to a mutable source BoxData with an optional shift.
*/
  template<class T, unsigned int C = 1, unsigned char D = 1, unsigned char E = 1>
  BoxData<T,C,D,E> alias(BoxData<T,C,D,E>& a_original, const Point& shift=Point::Zeros());

/// Alias (Const)
/**
   Creates a read-only alias to a const source BoxData with an optional shift.
*/
  template<class T, unsigned int C = 1, unsigned char D = 1, unsigned char E = 1>
  const BoxData<T,C,D,E> alias(const BoxData<T,C,D,E>& a_original, const Point& shift=Point::Zeros());

/// Slice Arbitrary Component (Non-Const)
/**
   Creates a BoxData<T,1,1,1> read-write alias to a prescribed component. 
*/
  template<typename T, unsigned int C, unsigned char D, unsigned char E>
  BoxData<T> slice(const BoxData<T,C,D,E>& a_src,
                   unsigned int a_c,
                   unsigned int a_d = 0,
                   unsigned int a_e = 0);

/// Slice Arbitrary Component Range (Non-Const)
/**
   Creates a BoxData<T,CC,1,1> read-write alias to a prescribed component range
*/
  template<typename T, unsigned int C, unsigned char CC>
  BoxData<T,CC,1,1> slice(const BoxData<T,C,1,1>& a_src, unsigned int a_nstart);


///@}
/** @name Pointwise Operators

   The suite of forall functions facilitate writing functions that operate pointwise on BoxData.
   To this end, the user must write a function with one of the following structures:
   @code
   PROTO_KERNEL_START
   void F_temp(Var<T,C,D,E>&, Args...)
   { ... }
   PROTO_KERNEL_END(F_temp, F)
   // OR
   PROTO_KERNEL_START
   void F_p_temp(Point&, Var<T,C,D,E>&, Args...)
   { ... }
   PROTO_KERNEL_END(F_p_temp, F_p)
   @endcode
   - PROTO_KERNEL_START and PROTO_KERNEL_END are required for cross-platform (e.g. GPU compatable) code.
   - The "#_temp" symbols are temporaries; the actual function symbol is the one without "_temp"
   - The template arguments of the first Var argument must match the output BoxData
   - The Point argument in the second signature corresponds to the Point of function application
   - Args... may include any number of Var& or read-only scalars.
   - The elements of Args... may have arbitrary tensor structure and const-ness
   - non-const objects in Args... have input-output semantics
   - The order and template arguments of the Vars in Args... must match the BoxData inputs of forall
   - If F is a member function of a class F MUST BE DECLARED STATIC
   - F or F_p may be an anonymous (lambda) function defined using the PROTO_LAMBDA macro

   Refer to the following code snippet for some sample valid forall input functions:
   \snippet Snippets.cpp proto_forall_func

*/
///@{

/// Pointwise Operator
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   This version of forall returns a new BoxData corresponding to the first argument of <code> a_F </code> which must be a <code> Var </code>.
   The domain of the created output is equal to the intersection of all other BoxData inputs (the function will fail if there are no other inputs).
   This function MUST specify the template parameters of the output BoxData. See the example code snippet below.

   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_constoprim
   Calling forall:
   \snippet Snippets.cpp proto_forall_1   

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
   \tparam T       (Mandatory!) Data type of return BoxData
   \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
   \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
   \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1

*/  
template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall(const Func& a_F, Srcs&&... a_srcs);

  ///same idea, but with flop counts and a timer name
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forallOp(unsigned long long int a_num_flops_point,   
                                   const char* a_timername, 
                                   const Func& a_F, Srcs&&... a_srcs);


/// Pointwise Operator: Overload with Box Argument
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   This version of forall returns a new BoxData corresponding to the first argument of <code> a_F </code> which must be a <code> Var </code>.
   This function MUST specify the template parameters of the output BoxData. See the example code snippet below.

   In general, this function should not be used unless absolutely necessary. Some valid use cases are:
   - Creating and initializing a BoxData without any BoxData inputs
   - Evaluating a pointwise function on a Box that is a proper subset of the intersection of all input BoxData

   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_constoprim
   Calling forall:
   \snippet Snippets.cpp proto_forall_2   

   \param a_F      Pointwise function
   \param a_box    Domain of computation. Must be a subset of the intersection of all input domains.
   \param a_srcs   Inputs (BoxData and primitive data)
   \tparam T       (Mandatory!) Data type of return BoxData
   \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
   \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
   \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall(const Func& a_F, Box  a_box, Srcs&&... a_srcs); 

///same idea, but with flop counts and a timer name
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forallOp(unsigned long long int a_num_flops_point,   
                                   const char* a_timername, 
                                   const Func& a_F, Box  a_box, Srcs&&... a_srcs); 

/// Pointwise Operator with Point Dependence
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   This version of forall allows the input function to be dependent on the Point at which it is applied.
   Hence, the first argument of a_F is a Point&, followed by the Var corresponding to the output BoxData
   This function MUST specify the template parameters of the output BoxData. See the example code snippet below.

   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_pointwise
   Calling forall:
   \snippet Snippets.cpp proto_forall_3   

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
   \tparam T       (Mandatory!) Data type of return BoxData
   \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
   \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
   \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall_p(const Func& a_F, Srcs&&... a_srcs);


///same idea, but with flop counts and a timer name
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forallOp_p(unsigned long long int a_num_flops_point,   
                                     const char* a_timername, 
                                     const Func& a_F, Srcs&&... a_srcs);

/// Pointwise Operator with Point Dependence: Overload with const Box Argument
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   This version of forall allows the input function to be dependent on the Point at which it is applied.
   Hence, the first argument of a_F is a Point&, followed by the Var corresponding to the output BoxData
   This function MUST specify the template parameters of the output BoxData. See the example code snippet below.

   In general, this function should not be used unless absolutely necessary. Some valid use cases are:
   - Creating and initializing a BoxData with only spatial dependence(e.g. F = f(x,y,z,dx,dy,dx)
   - Evaluating a pointwise function on a Box that is a proper subset of the intersection of all input BoxData

   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_pointwise
   Calling forall:
   \snippet Snippets.cpp proto_forall_4   

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
   \param a_box    Domain of computation. Must be a subset of the intersection of all input domains.
   \tparam T       (Mandatory!) Data type of return BoxData
   \tparam C       (Optional) Size of first component axis of return BoxData. Defaults to 1
   \tparam D       (Optional) Size of second component axis of return BoxData. Defaults to 1
   \tparam E       (Optional) Size of third component axis of return BoxData. Defaults to 1
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall_p(const Func& a_F, Box a_box, Srcs&&... a_srcs); 

///same idea, but with flop counts and a timer name
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forallOp_p(unsigned long long int a_num_flops_point,   
                                     const char* a_timername, 
                                     const Func& a_F, Box a_box, Srcs&&... a_srcs); 

/// In-Place Pointwise Operator
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   The "InPlace" versions of forall execute on existing BoxData and do not produce a new array.
   The domain of the created output is equal to the intersection of all BoxData inputs.

   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_constoprim
   Calling forall:
   \snippet Snippets.cpp proto_forall_5   

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace(const Func& a_F, Srcs&&... a_srcs);


///same idea, but with flop counts and a timer name
  template<typename Func, typename... Srcs>
  inline void forallInPlaceOp(unsigned long long int a_num_flops_point,   
                              const char* a_timername, 
                              const Func& a_F, Srcs&&... a_srcs);

/// In-Place Pointwise Operator on Prescribed Box
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   The "InPlace" versions of forall execute on existing BoxData and do not produce a new array.
   <code> a_F </code> will be applied at all points of the input <code> a_box </code>.

   In general, this function should not be used unless you want to restrict the domain of a_F's application
   to be something smaller than the intersection of all the inputs.   
 
   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_constoprim
   Calling forall:
   \snippet Snippets.cpp proto_forall_6   

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
   \param a_box    Domain of computation. Must be a subset of the intersection of all input domains.
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace(const Func& a_F, Box a_box, Srcs&&... a_srcs); 

///same idea, but with flop counts and a timer name
  template<typename Func, typename... Srcs>
  inline void forallInPlaceOp(unsigned long long int a_num_flops_point,   
                              const char* a_timername, 
                              const Func& a_F, Box a_box, Srcs&&... a_srcs); 


/// In-Place Pointwise Operator with Point Dependence
/**   
   Computes the function <code> a_F </code> at each Point of this BoxData.
   This version of forall allows the input function to be dependent on the Point at which it is applied.
   Hence, the first argument of <code> a_F </code> is a <code> Point& </code>, followed by the normal <code> Var </code> inputs

   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_pointwise
   Calling forall:
   \snippet Snippets.cpp proto_forall_7 

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace_p(const Func& a_F, Srcs&&... a_srcs);

///same idea, but with flop counts and a timer name
  template<typename Func, typename... Srcs>
  inline void forallInPlaceOp_p(unsigned long long int a_num_flops_point,   
                                const char* a_timername, 
                                const Func& a_F, Srcs&&... a_srcs);

/// In-Place Pointwise Operator with Point Dependence and Prescribed Box
/**
   Computes the function <code> a_F </code> at each Point of this BoxData.
   This version of forall allows the input function to be dependent on the Point at which it is applied.
   Hence, the first argument of <code> a_F </code> is a <code> Point& </code>, followed by the normal <code> Var </code> inputs

   In general, this function should not be used unless you want to restrict the domain of a_F's application
   to be something smaller than the intersection of all the inputs.   
 
   Example Usage:
   Input Function:
   \snippet Snippets.cpp proto_forall_pointwise
   Calling forall:
   \snippet Snippets.cpp proto_forall_8 

   \param a_F      Pointwise function
   \param a_srcs   Inputs (BoxData and primitive data)
   \param a_box    Domain of computation. Must be a subset of the intersection of all input domains.
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace_p(const Func& a_F,  Box a_box, Srcs&&... a_srcs); 

///same idea, but with flop counts and a timer name
  template<typename Func, typename... Srcs>
  inline void forallInPlaceOp_p(unsigned long long int a_num_flops_point,   
                                const char* a_timername, 
                                const Func& a_F,  Box a_box, Srcs&&... a_srcs); 

  // End of BoxData Doxygen Module
  /*@}*/

//========================================================================
//CUDA STUFF
//========================================================================
#ifdef PROTO_CUDA
  template<typename Func, typename... Srcs>
  inline void protoForall(const Func& a_F, Box a_box, Srcs&&... a_srcs); 



  template<typename Func, typename... Srcs>
  inline void protoForallStream(protoStream_t& a_stream, const Func& a_F, Box a_box, Srcs&&... a_srcs); 

  template<typename Func, typename... Srcs>
  inline void protoForall_p(const Func& a_F, Box a_box, Srcs&&... a_srcs); 
#endif

#include "implem/Proto_BoxDataImplem.H"
} //end Proto namespace
#endif //end include guard
