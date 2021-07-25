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
<<<<<<< HEAD
//enum MemType {HOST,DEVICE};
=======
>>>>>>> master
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
    
    /// Default constructor.
    inline LevelBoxData();

    /// Constructor over the DBL a_dbl. Calls define.
    inline LevelBoxData(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

    /// Define a default-constructed DisjointBoxLayout.
    inline void define(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

    /// NUmber of boxes in DBL.
    inline unsigned int size() const;

    /// Sets all the values on LBD to zero.
    inline void setToZero();

    /// fills ghost cells with valid data
    inline void exchange();
  
    ///fills ONLY VALID DATA
    /** call exchange if you want ghost cells filled */
    inline void copyTo(LevelBoxData& a_dest) const;

    /// gets the DBL.
    inline DisjointBoxLayout getDBL(){return m_dbl;};

      inline Point getGhostVector() {return m_ghost;}
    // DataIterator functions.

    /// Indexing.
    inline BoxData<T,C>& operator[](const DataIndex& dit);
    /// Const Indexing.
    inline const BoxData<T,C>& operator[](const DataIndex& dit) const;
    /// DataIterator accessed directly from the BLD. Begin iteration. 
    inline DataIterator begin() const;

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
