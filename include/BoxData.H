#ifndef _PROTO_RECTMDARRAY_H_
#define _PROTO_RECTMDARRAY_H_
#ifdef PROTO_CUDA
#include <thrust/transform_reduce.h>
#include <thrust/transform.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#endif

#include "Bx.H"
#include "PowerItoI.H"
#include "Stencil.H"

#include <memory> //for shared_ptr
#include <iostream>
#include <iomanip> //for pretty printing
#include <limits> // for max / min functions
#include <cstdlib> // for size_t
#include <cmath> // floating point abs
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
        // OR
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

/// Print with IOStream
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
    
    template <typename TT, unsigned int CC, unsigned char DD, unsigned char EE>
    friend class BoxData;
   
    template<typename Func, typename... TT>
    friend void pencilFunc(const Func& F,
                                    int count,
                                    TT&... vars);

    template<typename Func, typename... TT>
    friend void pencilFunc_p(const Func& F,
                                    int count,
                                    Point a_pt,
                                    TT&... vars);

    template<typename Func, typename... Srcs>
    friend void forallInPlaceBase(
                         const Func&            a_F,
                         Bx                     a_box,
                         Srcs&...               a_srcs); 
    
    template<typename Func, typename... Srcs>
    friend void forallInPlaceBase_p(
                         const Func&            a_F,
                         Bx                     a_box,
                         Srcs&...               a_srcs); 
    /// Pointwise Accessor
    /**
        Access component (c,d,e) of the const BoxData<T,C,D,E> associated with *this.
    */
    CUDA_DECORATION
    inline T& operator()(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0)
    {
        #ifdef __CUDA_ARCH__
        int idx = threadIdx.x + blockIdx.x*blockDim.x;
        T*  varloc =  (m_ptrs[a_c + C*a_d + C*D*a_e]);
        return varloc[idx];
        #else
        return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        #endif
    }


    /// Pointwise Accessor
    /**
        Access component (c,d,e) of the const BoxData<T,C,D,E> associated with *this.
    */
    CUDA_DECORATION
    inline const T& operator()(unsigned int a_c, unsigned char a_d = 0, unsigned char a_e = 0) const 
    {
        #ifdef __CUDA_ARCH__
        int idx = threadIdx.x + blockIdx.x*blockDim.x;
        return (m_ptrs[a_c + C*a_d + C*D*a_e])[idx];
        #else
        return *(m_ptrs[a_c + C*a_d + C*D*a_e]);
        #endif
    }
   
    //TODO increment operations should be private with necessary friend additions -clg 
    inline Var& operator++() 
    {
        for (int ii = 0; ii < C*D*E; ii++)
        {
            ++m_ptrs[ii];
        }
        return *this;
    }
    inline const Var& operator++() const
    {
        for (int ii = 0; ii < C*D*E; ii++)
        {
            ++m_ptrs[ii];
        }
        return *this;
    }
    inline Var& operator--() 
    {
        for (int ii = 0; ii < C*D*E; ii++)
        {
            --m_ptrs[ii];
        }
        return *this;
    }
    
    inline const Var& operator--() const
    {
        for (int ii = 0; ii < C*D*E; ii++)
        {
            --m_ptrs[ii];
        }
        return *this;
    }

private:
   T* m_ptrs[C*D*E];
};



/// Multidimensional Rectangular Array
/**
    BoxData is the main dataholder class of Proto. It contains an 
    array of data of type T on a domain defined by a Bx.
    Each data point is associated with a Point in the Bx.
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
  
    /// Default Constructor
    /**
       Allocates no space; user must call define on resulting object.
    */
    BoxData();

    /// Box Constructor
    /**
      Creates and allocates an uninitialized array based on a Bx input. 
      \param a_box  Bx defining the domain of *this
    */
    BoxData(const Bx& a_box);
    
    /// Box Constructor With Initialization
    /**
       Identical to calling the Box constructor followed by setVal(a_init).
       Not recommended in performance-critical code if initialization is unnecessary.
        
       \param a_box    Bx defining the domain of *this
       \param a_init   T value used to initialize all data values
    */
    BoxData(const Bx& a_box, T a_init);
    
    /// Copy constructor
    /**
       Performs a deep copy of all data in src
       
       \param a_src   Source data
    */
    BoxData(const BoxData<T,C,D,E>& a_src);

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
       
       \param a_box   Bx defining the domain of *this
    */  
    void define(const Bx& a_box);
   
    //TODO: It does not make sense to use a_ncomp here. We should try to get rid of that -clg 
    /// Raw Pointer Constructor
    /**
       Builds a vector valued BoxData by aliasing a raw pointer to T which should be of size
       <code> a_box.size()*a_ncomp </code>. The data in <code> a_ptr </code> is not copied.
       This constructor is not recommended for public use, but is provided to facilitate
       compatability with third party libraries.
       
       \param a_ptr   Source buffer of type T to which *this will be aliased
       \param a_box   Bx defining the domain of *this
       \param a_ncomp (Optional) Dummy input for the number of components. Completely unused.
    */
    BoxData(const T* a_ptr, const Bx& a_box, int a_ncomp = C);

    //TODO: It does not make sense to use a_ncomp here. We should try to get rid of that -clg 
    /// Raw Pointer Define
    /**
       Defines a vector valued BoxData by aliasing a raw pointer to T which should be of size
       <code> a_box.size()*a_ncomp </code>. The data in <code> a_ptr </code> is not copied.
       This constructor is not recommended for public use, but is provided to facilitate
       compatability with third party libraries.
       
       \param a_ptr   Source buffer of type T to which *this will be aliased
       \param a_box   Bx defining the domain of *this
       \param a_ncomp (Optional) Dummy input for the number of components. Completely unused.
    */
    void define(const T* a_ptr, const Bx& a_box, int a_ncomp = C);

    //TODO: Not really part of the public API. This should maybe be private -clg
    /// LazyStencil Constructor
    /**
        Builds a BoxData from a lazily evaluated Stencil operation. Called implicitly.
    
        \param a_op   A LazyStencil created by Proto::Stencil::operator() or Proto::InterpStencil::operator()
    */
    BoxData(LazyStencil<T,C,D,E>&& a_op);
   
    //TODO: Should be private? 
    /// Slice Constructor
    /**
        Used internally. Not part of the public interface.
    */
    BoxData(shared_ptr<T> a_data,const T* a_ptr,const Bx& a_box);


    /// Destructor.
    ~BoxData(); 
  
    ///@}
   
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Data Movement */
    ///@{

    /// Copy Assignment Operator
    /**
       Performs a deep copy

       \param a_src   Source data
    */
    BoxData& operator=(const BoxData<T,C,D,E>& a_src);

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
       Copy with a prescribed Bx argument and optional shift. Explicitly, this function copies 
       the subset of data from <code> *this </code> contained in <code> a_srcBox </code>
       into the region of <code> a_dest </code> defined by <code> a_srcBox.shift(a_destShift) </code>.

       This function will fail if <code> a_srcBox </code> is not contained in both <code> this->box() </code>
       and <code> a_dest.box().shift(a_destShift) </code>.

       \param a_dest       Destination data holder
       \param a_srcBox     Region of data to copy from *this
       \param a_destShift  (Optional) Determines region of a_dest to copy data to.
    */
    void copyTo(BoxData<T,C,D,E>&   a_dest,
                const Bx&           a_srcBx,
                const Point&        a_destShift = Point::Zeros()) const;
    

    /// General Copy                                                                                                                                                                                            
    /**
        The most general form of copyTo. Copies a prescribed set of components from <code> *this </code>
        into <code> a_dest </code> within prescribed region with a possible shift. 
       
        Example Usage:
        \snippet Snippets.cpp proto_copyTo 
  
        \param a_dest       Destination data holder
        \param a_srcBox     Region of data to copy from *this
        \param a_srcComps   Components of *this to copy from                                                            
        \param a_destShift  Determines region of a_dest to copy data to (e.g. a_srcBox.shift(a_destShift))
        \param a_destComps  Components of a_dest to copy into. Must be the same size as a_srcComps.                                                                                                             
    */
    template <unsigned int Cdest,
              unsigned char Ddest,
              unsigned char Edest>
    void copyTo(BoxData<T,Cdest,Ddest,Edest>& a_dest,
                const Bx&    a_srcBox,
                CInterval    a_srcComps,
                const Point& a_destShift,
                CInterval    a_destComps) const;

    // TODO: This ought to live in Chombo in my opinion... -clg
    /// General Copy From
    /**
        Nearly identical to the most general version of <code> copyTo </code> except the copy direction is reversed.
        Only valid for vector valued BoxData. This function is provided to facilitate interoperability with third party libraries.
    
        \param a_dsrc     Source data
        \param a_srcBox   Source domain to copy from
        \param a_srcComp  First component to copy from source
        \param a_destBox  Destination domain to copy into
        \param a_destComp First component to copy into destiation
        \param a_numcomp  Number of components to copy
    */
    template<unsigned int Csrc>
    void copy(const BoxData<T,Csrc,D,E>& a_dsrc,
              const Bx&                  a_srcBox,
              unsigned int               a_srcComp,
              const Bx&                  a_destBox,
              unsigned int               a_destComp,
              unsigned int               a_numcomp);
    
    ///@}
  
    ///////////////////////////////////////////////////////////////////////////////////////////////
    /** @name Accessors */
    ///@{
    
    /// Read Only Data Accessor 
    /**
      Read-only access to data stored in <code> *this </code>. This function is read-only to facilitate
      cross-platform code (e.g. GPU compatability). Pointwise write operations should be done through forall.

      \param a_pt   A Point in the domain of <code> *this </code>
      \param a_c    Value of first tensor index
      \param a_d    Value of second tensor index
      \param a_e    Value of third tensor index
    */
    CUDA_DECORATION
    inline T operator()(const Point& a_pt,
                        unsigned int  a_c = 0,
                        unsigned char a_d = 0,
                        unsigned char a_e = 0) const
    {
      PROTO_ASSERT(m_box.contains(a_pt),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_pt not in this->box()");
      PROTO_ASSERT((a_c >= 0) && (a_c < C),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_c = %i. a_c must be in [0,%i)",a_c,C);
      PROTO_ASSERT((a_d >= 0) && (a_d < D),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_d = %i. a_d must be in [0,%i)",a_d,D);
      PROTO_ASSERT((a_e >= 0) && (a_e < E),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e) invalid for a_e = %i. a_e must be in [0,%i)",a_e,E);
      return m_rawPtr[index(a_pt,a_c,a_d,a_e)];
    }
  
    //TODO: Should be private? 
    /// Index Accessor (Non-Const)
    /**
       Access *this using a linear index. Useful for interacting with
       BoxData like a regular buffer, but not recommended.

       \param a_index   Index in [0,this->size() = m_box.size()*C*D*E)
    */
    inline T* operator[](unsigned int a_index); 
    
    //TODO: Should be private? 
    /// Index Accessor (Const)
    /**
        Access *this using a linear index. Useful for interacting with
        BoxData like a regular buffer, but not recommended.

        \param a_index   Index in [0,this->size() = m_box.size()*C*D*E)
    */
    inline const T* operator[](unsigned int a_index) const; 
   
    //TODO: Should be private? 
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
      unsigned int ioff = index(a_p, a_c, a_d, a_e);
      return &(m_rawPtr[ioff]);
    }

    //TODO: Should be private? 
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
      unsigned int ioff = index(a_p, a_c, a_d, a_e);
      return &(m_rawPtr[ioff]);
    }
    
    //TODO: Should be private? -clg
    /// Access Shared Pointer
    /**
        Used internally. Not recommended for public use. 
    */  
    inline shared_ptr<T> getData() const
    {
      return m_data; ///< Data array
    }
   
    
    //used by slicing
    /* TODO: Flagged for removal; not used anymore I don't think -clg
    T* getRawPtr()
    {
      return m_rawPtr;
    } 
    const T* getRawPtr() const 
    {
      return m_rawPtr;
    } 
    */

    /// Compute Index
    /**
        Computes the index in [0,this->size()) associated with a given Point and components
      
        \param a_pt   A Point in the domain of <code> *this </code>
        \param a_c    Value of first tensor index
        \param a_d    Value of second tensor index
        \param a_e    Value of third tensor index
    */
    CUDA_DECORATION
    inline unsigned int index(const Point a_pt,
                              unsigned int a_c = 0,
                              unsigned int a_d = 0,
                              unsigned int a_e = 0) const
    {
      PROTO_ASSERT(m_box.contains(a_pt),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
    invalid for a_pt not in this->box()");
      PROTO_ASSERT((a_c >= 0) && (a_c < C),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
    invalid for a_c = %i. a_c must be in [0,%i)",a_c,C);
      PROTO_ASSERT((a_d >= 0) && (a_d < D),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
    invalid for a_d = %i. a_d must be in [0,%i)",a_d,D);
      PROTO_ASSERT((a_e >= 0) && (a_e < E),
                   "BoxData::operator(Point a_pt, uint a_c, uchar a_d, uchar a_e)\
    invalid for a_e = %i. a_e must be in [0,%i)",a_e,E);
      int m = m_box.size();
      int k = m_box.index(a_pt);
      return k+m*a_c+a_d*m*C+a_e*m*C*D;
    }
  
    // Get Domain Box
    /**
       Return's the domain Bx of <code> *this </code>. 
       Includes ghost cells if it was built it that way.
    */
    inline Bx box() const {return m_box;};
    
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
   
    //TODO: Should definitely not be public 
    /// Create Pointwise Access Variable (Non-Const)
    /**
       Not for public use
    */
    inline Var<T,C,D,E> var(const Point& a_pt);
    
    //TODO: Should definitely not be public 
    /// Create Pointwise Access Variable (Const)
    /**
       Not for public use
    */
    inline Var<T,C,D,E> var(const Point& a_pt) const;
 
#ifdef PROTO_CUDA
    //====================================================================
    //  CUDA HELPER FUNCTIONS
    //====================================================================
    // TODO: flagged for removal -clg
    //inline Var<T,C,D,E> cudaVar(const Point& a_pt);
    
    // TODO: flagged for removal -clg
    //inline Var<T,C,D,E> cudaVar(const Point& a_pt) const;

#endif
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
                const Bx& a_box);

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
                const Bx& a_box,
                int a_c,
                int a_d = 0,
                int a_e = 0);

    /// Absolute Maximum Value (Global)
    /**
       Returns the maximum absolute value over the entire data set
    */
    T absMax() const;

    /// Absolute Maximum Value (Componentwise)
    /**
       Returns the maximum absolute value of a given component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T absMax(int a_c, int a_d = 0, int a_e = 0) const;

    /// Minimum Value (Global)
    /**
       Returns the minimum value over the entire data set
    */
    T min() const;

    /// Minimum Value (Componentwise)
    /**
       Returns the minimum value of a given component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T min(int a_c, int a_d = 0, int a_e = 0) const;

    /// Maximum Value (Global)
    /**
       Returns the maximum value over the entire data set
    */
    T max() const;

    /// Maximum Value (Componentwise)
    /**
       Returns the maximum value of a given component
       
       \param a_c      First tensor index.
       \param a_d      Second tensor index.
       \param a_e      Third tensor index.
    */
    T max(int a_c, int a_d = 0, int a_e = 0) const;

    /// Shift Domain
    /**
       Shifts the domain Bx of this, moving all data along with it.
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
                   const Bx& a_box,
                   CInterval a_comps) const;

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
                  const Bx& a_box,
                  CInterval a_comps);

   
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
    void printData(const Bx& a_box, int a_prec = 2) const;

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
    void printData(const Bx& a_box, int a_c, int a_d, int a_e, int a_prec = 2) const;

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
                                            const Point&                    a_shift);

    template<typename TT, unsigned int CC, unsigned char DD, unsigned char EE>
    friend BoxData<TT> slice(BoxData<TT,CC,DD,EE>& a_src,
                             unsigned int a_c,
                             unsigned int a_d,
                             unsigned int a_e);
  private:
    
    Bx                      m_box; ///< Bx defining the domain of *this
    ::std::shared_ptr<T>    m_data; ///< Data array
    T*                      m_rawPtr; ///< Raw pointer to the data
  };

///////////////////////////////////////////////////////////////////////////////////////////////
/** @name Alias and Slice Operators */
///@{

/// Alias (Non-Const)
/**
   Create an alias to a source (mutable) BoxData with an optional shift.

   Example:
   @code
   //DIM=2
   using namespace Proto;
   Box B = Box::Cube(4); //[(0,0), (3,3)]
   Point shift = Point::Basis(0);
   BoxData<double> R0(B,1337);
   auto R1 = alias(R0,shift);
   R0.box().print(); //prints [(0,0), (3,3)]
   R1.box().print(); //prints [(1,0), (4,3)]
   for (auto iter = B.begin(); iter != B.end(); ++iter)
   {
   std::cout << (&R0(*iter) == &R1(*iter + shift)) << std::endl; //true for all *iter
   }
   // or alternatively:
   for (int ii = 0; ii < R0.size(); ii++)
   {
   std::cout << (&R0[ii] == &R1[ii]) << std::endl; //true for all ii
   }
   @endcode

*/
  template<class T, unsigned int C = 1, unsigned char D = 1, unsigned char E = 1>
  BoxData<T,C,D,E> alias(BoxData<T,C,D,E>& a_original, const Point& shift=Point::Zeros());

/// Alias (Const)
/**
   Create an alias to a source (immutable) BoxData with an optional shift.
   See documentation of the non-const version for example usage.
*/
  template<class T, unsigned int C = 1, unsigned char D = 1, unsigned char E = 1>
  const BoxData<T,C,D,E> alias(const BoxData<T,C,D,E>& a_original, const Point& shift=Point::Zeros());

/// Slice Arbitrary Component (Non-Const)
/**
   Creates a BoxData<T,1,1,1> alias to a prescribed component. 

   Example:
   @code
   //DIM=2
   using namespace Proto;
   Box B = Box::Cube(4); //[(0,0), (3,3)]
   BoxData<double,3,3> R0(B);
   // R1 is a slice of only the (2,1) component
   auto R1 = slice(R0,2,1);
   for (auto iter = B.begin(); iter != B.end(); ++iter)
   {
   std::cout << (&R0(*iter,2,1) == &R1(*iter)) << std::endl; //prints true for all *iter
   }
   @endcode

*/
  template<typename T, unsigned int C, unsigned char D, unsigned char E>
  BoxData<T> slice(BoxData<T,C,D,E>& a_src,
                   unsigned int a_c,
                   unsigned int a_d = 0,
                   unsigned int a_e = 0);

/// Slice Arbitrary Component (Const)
/**
   Creates a BoxData<T,1,1,1> alias to a prescribed component.
   See the documentation for the non-const version for example usage. 
*/
  template<typename T, unsigned int C, unsigned char D, unsigned char E>
  const BoxData<T> slice(const BoxData<T,C,D,E>& a_src,
                         unsigned int a_c,
                         unsigned int a_d = 0,
                         unsigned int a_e = 0);


/** @name Pointwise Operators*/
///@{

/// Pointwise Operator
/**
   Evaluate a pointwise function Func.
   Func must have one of the following signatures:

   @code
   void F(Var<T,C,D,E>&, Args...);
   // OR
   void F_p(Point&, Var<T,C,D,E>&, Args...); 
   @endcode
   - The template arguments of the first Var argument must match the output BoxData
   - The Point argument in the second signature corresponds to the Point of function application
   - Args... may include any number of Var& or Param& objects in any order
   - The template arguments of the Vars and Params in Args... are arbitrary
   - Objects in Args... may be const or non-const as needed
   - non-const objects in Args... have input-output semantics
   - The order and template arguments of the Vars in Args... must match the BoxData inputs
   - If F is a member function of a class F MUST BE DECLARED STATIC

   If a Bx parameter is provided, the returned BoxData will be defined based on it. Otherwise,
   the returned BoxData will be defined on the intersection of all input BoxData. In the latter case,
   the function will fail if no input BoxData are specified.

   \param a_F      Pointwise function
   \param a_srcs   Parameters
   \tparam T       Data type of return BoxData
   \tparam C       Size of first component axis of return BoxData. Defaults to 1
   \tparam D       Size of second component axis of return BoxData. Defaults to 1
   \tparam E       Size of third component axis of return BoxData. Defaults to 1

   Usage:
   @code
    
   using namespace Proto;

   // input function for forall
   void foo(Var<double>& KE, const Var<double, DIM>& U, Param<double>& KE_max)
   {
   double u_squared = 0.0;
   for (int dir = 0; dir < DIM; dir++)
   {
   u_squared += U(dir)*U(dir); //access BoxData in component dir through U
   }
   double ke = u_squared*0.5;
   KE_max = (ke > KE_max) ? ke : KE_max; //access Param value KE_max like normal POD
   KE(0) = ke;
   }

   int main(int argc, char** argv)
   {
   Bx b = Bx::Cube(16);
   BoxData<double, DIM> U(b);
   double ke_max = 0.0;
   // ... initialize U ... 
   // call forall. Arguments must be in the same order as in foo.
   // The first argument from foo - KE - corresponds to the return value
   // No copy occurs here due to C++ move semantics
   // KE will be defined using the Bx of U, namely b. 
   BoxData<double> KE = forall<double>(foo, U, ke_max)
   }
   @endcode

   All versions of forall are also compatible with lambda (anonymous) functions, as shown below:

   @code
   int main(int argc, char** argv)
   {
   Bx b = Bx::Cube(16);
   BoxData<double, DIM> U(b);
   double ke_max = 0.0;
   // ... initialize U ... 

   // Call to forall with lambda function
   BoxData<double> KE = forall<double>(
   [](Var<double>& a_ke, const Var<double, DIM>& a_u, Param<double> a_ke_max) //lambda start
   {
   double u_squared = 0.0;
   for (int dir = 0; dir < DIM; dir++)
   {
   u_squared += a_u(dir)*a_u(dir);
   }
   double ke = u_squared*0.5;
   a_ke_max = (ke > KE_max) ? ke : a_ke_max; 
   a_ke(0) = ke;
   },                                                                          //lambda end
   U, ke_max);                                                                 //remaining arguments
   }
   @endcode
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall(const Func& a_F, Srcs&... a_srcs);

/// Pointwise Operator: Overload with Box Argument
/**
   Supplying an invalid Bx argument may result in invalid code.
   this function is only recommended when it is absolutely necessary.
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall(const Func& a_F, Bx  a_box, Srcs&... a_srcs); 

/// Pointwise Operator with Point Dependence
/**
   All versions of forall with _p may be called with a function using the signature:
    
   @code
   void F_p(Point&, Var<T,C,D,E>&, Args...); 
   @endcode

   The additional Point argument allows the user to do pointwise computations that are
   dependent on the point of computation. For example:

   @code
   using namespace Proto;
   void square(Point& a_p, Var<double>& a_y)
   {
   double x = (double)a_p[0]; //access the current Point;
   a_y(0) = x*x;
   }

   int main(int argc, char** argv)
   {
   int domainSize = 16;
   Bx b = Bx::Cube(domainSize);
   BoxData<double> X(b);
   //initialize X;
   auto Y = forall_p<double>(square, X);
   }

*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall_p(const Func& a_F, Srcs&... a_srcs);

/// Pointwise Operator with Point Dependence: Overload with const Box Argument
/**
   Supplying an invalid Bx argument may result in invalid code.
   this function is only recommended when it is absolutely necessary.
*/
  template<typename T, unsigned int C=1, unsigned char D=1, unsigned char E=1,
           typename Func, typename... Srcs>
  inline BoxData<T,C,D,E> forall_p(const Func& a_F, Bx a_box, Srcs&... a_srcs); 

/// In-Place Pointwise Operator (Old Forall)
/**
   Modified version of forall which does not return a BoxData.
   All arguments have input-output semantics.
    
   @code
   using namespace Proto;

   // input function for forall
   void foo(Var<double>& KE, const Var<double, DIM>& U, Param<double>& KE_max)
   {
   double u_squared = 0.0;
   for (int dir = 0; dir < DIM; dir++)
   {
   u_squared += U(dir)*U(dir); //access BoxData in component dir through U
   }
   double ke = u_squared*0.5;
   KE_max = (ke > KE_max) ? ke : KE_max; //access Param value KE_max like normal POD
   KE(0) = ke;
   }

   int main(int argc, char** argv)
   {
   Bx b = Bx::Cube(16);
   Bx B = b.grow(1); //B is slightly larger than b;
   BoxData<double, DIM> U(b);
   BoxData<double> KE(B,0);
   double ke_max = 0.0;
   // template arguments are unnecessary for in-place forall
   forallInPlace(foo, KE, U, ke_max);
   // because KE and U have different Bx domains (B and b) forall will compute foo
   // on the intersection (namely b).
   }
   @endcode
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace(const Func& a_F, Srcs&... a_srcs);

/// In-Place Pointwise Operator (Old Forall): Overload with const Box Argument
/**
   Supplying an invalid Bx argument may result in invalid code.
   this function is only recommended when it is absolutely necessary.
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace(const Func& a_F, Bx a_box, Srcs&... a_srcs); 

/// In-Place Pointwise Operator with Point Dependence
  template<typename Func, typename... Srcs>
  inline void forallInPlace_p(const Func& a_F, Srcs&... a_srcs);

/// In-Place Pointwise Operator with Point Dependence: Overload with Box Argument
/**
   Supplying an invalid Bx argument may result in invalid code.
   this function is only recommended when it is absolutely necessary.
*/
  template<typename Func, typename... Srcs>
  inline void forallInPlace_p(const Func& a_F,  Bx a_box, Srcs&... a_srcs); 

//========================================================================
//CUDA STUFF
//========================================================================
#ifdef PROTO_CUDA
  template<typename Func, typename... Srcs>
  inline void cudaForall(const Func& a_F, Bx a_box, Srcs&... a_srcs); 

  template<typename Func, typename... Srcs>
  inline void cudaForall_p(const Func& a_F, Bx a_box, Srcs&... a_srcs); 
#endif

#include "BoxDataImplem.H"
} //end Proto namespace
#endif //end include guard
