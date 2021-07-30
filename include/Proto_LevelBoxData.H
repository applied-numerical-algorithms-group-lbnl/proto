#ifndef _PROTO_LEVELBOXDATA_H_
#define _PROTO_LEVELBOXDATA_H_

#include "Proto_MemType.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_Copier.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelIterators.H"
#include "Proto_WriteBoxData.H"
#include "Proto_MayDay.H"
#include "Proto_SPMD.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cstring> // for Writing data to .vtk files.
#define LBD LevelBoxData<T,C,MEMTYPE,CENTERING>
#define LDOP LDOperator<T,C,MEMTYPE,CENTERING>
using namespace std;
using std::shared_ptr;

namespace Proto 
{
  template
  <typename T, unsigned int C, MemType MEMTYPE, unsigned int CENTERING >
  class LevelBoxData;
    template <typename T,unsigned int C=1, MemType MEMTYPE = MemType::HOST ,unsigned int CENTERING=0 >
  class LDOperator
  {
  public:
    LBD* m_srcPtr;
    LBD* m_destPtr;
    ~LDOperator()
    {
    }
    inline LDOperator(LBD& a_src,LBD& a_dest)
    {
      m_srcPtr = &a_src;
      m_destPtr = &a_dest;
    }
    inline int size(const Box& a_bx) const
    {
      return C*a_bx.size()*sizeof(T);
    }
    inline void linearOut(void* a_buf, 
                          const Box& a_bx,
                          const DataIndex& a_diFr
                           ) const
    {
      BoxData<T,C>& bd = m_srcPtr->operator[](a_diFr);
      CInterval cint(0,C-1);
      bd.linearOut(a_buf, a_bx,cint);
    }
    inline void linearIn(void* a_buf, 
                         const Box& a_bx,
                         const DataIndex& a_diTo) const
    {
      BoxData<T,C>& bd = m_destPtr->operator[](a_diTo);
      CInterval cint(0,C-1);
      bd.linearIn(a_buf, a_bx,cint);
    }

    inline void op(const Box& a_bxFr,const DataIndex& a_diFr,
                   const Box& a_bxTo,const DataIndex& a_diTo) const
    {
      int LinearizationTest = 1;
      if (LinearizationTest == 0)
        {
          PR_TIMERS("LOP_localCopy_1");
          BoxData<T,C>& src = m_srcPtr->operator[](a_diFr);
          BoxData<T,C>& dest = m_destPtr->operator[](a_diTo);
          Point shift = a_bxTo.low() - a_bxFr.low();
          
          src.copyTo(dest, a_bxFr, shift);
        }
      else
        {
           PR_TIMERS("LOP_localCopy_2");
          int sizeSource = size(a_bxFr);
          int sizeDest   = size(a_bxTo);
          
          if (sizeSource != sizeDest)
            {
              MayDay<void>::Abort("LinearizationTest failure: dest and source have different sizes");
            }
          
          vector<char> buffer(sizeSource);
          {
            PR_TIMERS("linearIn / linearOut inside LDOP::op");
            void* b = (void*)&(buffer[0]);
            linearOut(b, a_bxFr,a_diFr);
            linearIn(b, a_bxTo,a_diTo);
          }
        }
    }
  };
///data over a DisjointBoxLayout with ghost cells.
  template
  <typename T, unsigned int C=1, MemType MEMTYPE=MemType::HOST, unsigned int CENTERING=DIM >
  class LevelBoxData
  {
  public:
    
    /// Default Constructor
    inline LevelBoxData();

    /// Layout Constructor
    /**
        Builds a non-trivial LevelBoxData from a DisjointBoxLayout with specified ghost cells
        Calls LevelBoxData::define(...)
    */
    inline LevelBoxData(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

    /// Lazy Constructor
    /**
        Lazily construct or redefine a LevelBoxData
        FIXME: does this function handle previously defined LBDs?
    */
    inline void define(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

    /// Size (in Patches)
    /**
        Returns the number of patches in *this
        FIXME: is this function local or global?
        FIXME: rename this function so it doesn't collide with the serial size. (numPatches? numBoxes?)
    */
    inline unsigned int size() const;

    /// Sets all the values on BLD to zero.
    inline void setToZero();

    /// fills ghost cells with valid data
    inline void exchange();
  
    ///fills ONLY VALID DATA
    /** call exchange if you want ghost cells filled */
    inline void copyTo(LevelBoxData& a_dest) const;
    
    /// Linear In (Patch, Box)
    /**
        Reads data from a serial buffer, populating the patch associated with 
        a_index within the subset specified by a_box. 

        \param a_buffer     The read buffer
        \param a_box        A box defining the subset to read into
        \param a_index      The index defining the patch to read into
    */ 
    inline void linearIn(void* a_buffer, 
                         const Box& a_box,
                         const DataIndex& a_index);
    
    /// Linear Out (Patch, Box)
    /**
        Writes data to a serial buffer, from the patch associated with 
        a_index within the subset specified by a_box. 
        
        \param a_buffer     The write buffer
        \param a_box        A box defining the subset to write from
        \param a_index      The index defining the patch to write from
    */ 
    inline void linearOut(void* a_buffer, 
                         const Box& a_box,
                         const DataIndex& a_index) const;
    
    /// Serial Size (Patch, Box)
    /**
        Computes the size in bytes of the serial data buffer needed to hold
        the patch associated with a_index within the subset specified by a_box.
        
        \param a_box        A box defining the subset of a patch
        \param a_index      The index defining the patch
    */ 
    inline size_t size(const Box& a_box,
                       const DataIndex& a_index) const;
    /// Linear In (Patch)
    /**
        Overload of LinearIn that always reads in a full patch
    */
    inline void linearIn(void* a_buffer, 
                         const DataIndex& a_index);
    /// Linear Out (Patch)
    /**
        Overload of LinearOut that always writes out a full patch
    */
    inline void linearOut(void* a_buffer, 
                         const DataIndex& a_index) const;
    /// Size (Patch)
    /**
        Overload of Size that always computes the size of a full patch
    */
    inline size_t size( const DataIndex& a_index) const;
    
    /// Linear In (All Local Data)
    /**
        Overload of LinearIn that reads all local data
    */
    inline void linearIn(void* a_buffer);
    /// Linear Out (All Local Data)
    /**
        Overload of LinearOut that writes out all local data
    */
    inline void linearOut(void* a_buffer) const;
    /// Size (All Local Data)
    /**
        Overload of Size that computes the full local size of the data in *this
        FIXME: This function collides with the "number of boxes" version of size. 
    */
    inline size_t serialSize() const;

    /// gets the DBL.
    inline DisjointBoxLayout layout(){return m_dbl;}; // more expressive -clg
    inline DisjointBoxLayout getDBL(){return m_dbl;};

    inline Point ghost(){return m_ghost;} //more expressive -clg
    inline Point getGhostVector() {return m_ghost;}
    // DataIterator functions.

    /// Indexing.
    inline BoxData<T,C>& operator[](const DataIndex& dit);
    /// Const Indexing.
    inline const BoxData<T,C>& operator[](const DataIndex& dit) const;
    /// DataIterator accessed directly from the BLD. Begin iteration. 
    inline DataIterator begin() const;

    /// Write to HDF5
    template<typename... Args>
    inline void write(array<double, DIM> a_dx,
        std::string a_filename,
        Args... a_params) const;
    
    /// Write to HDF5
    template<typename... Args>
    inline void write(std::vector<std::string> a_varNames,
        array<double, DIM> a_dx,
        std::string a_filename,
        Args... a_params) const;
    
    /// Write to HDF5
    template<typename... Args>
    inline void write(std::vector<std::string> a_varNames,
        std::string a_filename,
        Args... a_params) const;
    
    /// Write to HDF5
    template<typename... Args>
    inline void write(std::string a_filename,
        Args... a_params) const;

    /// For debugging purposes.
    int s_verbosity = 0;
  private: 

    std::vector<shared_ptr<BoxData<T,C> > > m_data;
    Point                                   m_ghost;
    DisjointBoxLayout                       m_dbl;
    bool                                    m_isDefined;
    Copier                                  m_exchangeCopier;
    // disallow copy constructors and assignment operators
    // to avoid very hard-to-find performance problems
    LBD& operator= (const LBD& a_rhs);
    LBD(const LBD& a_rhs);
  /// Stand-alone function that writes a LevelData on a domain-decomposed rectangle into a .vtk file.
  // inline 
  //void WriteData(LevelBoxData<BoxData<double>& a_phi,int a_iter,double a_dx,std::string a_str);

protected:
  // Internal fucntions that execute copying, exchange.
 void makeItSo(
               const Copier&       a_copier,
               const LDOP&         a_op) const;
  
  void makeItSoBegin(
                     const Copier&       a_copier,
                     const LDOP&         a_op) const; 

  void makeItSoLocalCopy(
                         const Copier&       a_copier,
                         const LDOP&         a_op) const;

  void makeItSoEnd(const LDOP&         a_op) const;
  
  //========================================================================
  //
  // data structures used by makeItSo when we have some
  // data that needs to be moved (ie. there are entries
  // in the 'FROM' or 'TO' CopyIterators)
  //
  void completePendingSends() const;

  void allocateBuffers(
                       const Copier&   a_copier,
                       const LDOP& a_op) const;

  void writeSendDataFromMeIntoBuffers(const LDOP& a_op) const;
  void postSendsFromMe() const ;

  void postReceivesToMe() const ;

  void unpackReceivesToMe(const LDOP& a_op) const;

  mutable CopierBuffer* m_buff;
  };
#include "implem/Proto_LevelBoxDataImplem.H"
} //end namespace Proto
#endif
