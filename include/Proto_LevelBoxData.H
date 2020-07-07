#ifndef _PROTO_LEVELBOXDATA_H_
#define _PROTO_LEVELBOXDATA_H_

#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_Copier.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelIterators.H"
#include "Proto_WriteBoxData.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cstring> // for Writing data to .vtk files.
enum MemType {HOST,DEVICE};
#define LBD LevelBoxData<T,C,MEMTYPE,CENTERING>
#define LDOP LDOperator<T,C>
using namespace std;
using std::shared_ptr;

namespace Proto 
{
  template <typename T,unsigned int C=1>
  class LDOperator
  {
  public:
    virtual ~LDOperator()
    {
    }
    virtual int size(BoxData<T,C>& a_bd, const Box& a_bx) const
    {
      //    return arg.size(b, comps);
      return a_bd.charsize(a_bx);
    }
    virtual void linearOut(BoxData<T,C>& a_bd, 
                           void* a_buf, 
                           const Box& a_bx
                           ) const
    {
      a_bd.linearOut(a_buf, a_bx);
    }
    virtual void linearIn(BoxData<T,C>& a_bd,  
                          void* a_buf, 
                          const Box& a_bx) const
    {
      a_bd.linearIn(a_buf, a_bx);
    }

    virtual void op(BoxData<T,C>& a_dest,
                    const Box& a_bxFr,
                    const Box& a_bxTo,
                    const BoxData<T,C>& a_src) const
    {
      if (LinearizationTest == 0)
        {          
          a_dest.copyTo(a_src, a_bxFr, a_bxTo);
        }
      else
        {
          int sizeSource = size(a_src,  a_bxFr);
          int sizeDest   = size(a_dest, a_bxTo);

          
            
          if (sizeSource != sizeDest)
            {
              MayDay::Abort("LinearizationTest failure: dest and source have different sizes");
            }
            
          vector<char> buffer(sizeSource);
          void* b = (void*)&(buffer[0]);
          linearOut(a_src, b, a_bxFr);
          linearIn(a_dest, b, a_bxTo);
        }
    }
  };
///data over a disjointboxlayout with ghost cells.
  template
  <typename T, unsigned int C=1, MemType MEMTYPE=HOST, unsigned int CENTERING=DIM >
  class LevelBoxData
  {
  public:
    
    ///
    inline LevelBoxData();

    /// 
    inline LevelBoxData(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

    /// 
    inline void define(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

    ///
    inline unsigned int size() const;

    ///
    inline void setToZero();

    /// fills ghost cells with valid data
    inline void exchange();
  
    ///fills ONLY VALID DATA
    /** call exchange if you want ghost cells filled */
    inline void copyTo(LevelBoxData& a_dest) const;

    /// gets the dbl PC: added.
    inline DisjointBoxLayout getDBL(){return m_dbl;};
    // DataIterator functions.

    /// Indexing.
    inline BoxData<T,C>& operator[](const DataIndex& dit);
    /// Const Indexing.
    inline const BoxData<T,C>& operator[](const DataIndex& dit) const;
    /// Begin iteration. 
    inline DataIterator begin() const;

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
 void makeItSo(const LBD&          a_src,
               LBD&                a_dest,
               const Copier&       a_copier,
               const LDOP&         a_op =LDOP()) const;
  
  void makeItSoBegin(const LBD&          a_src,
                     LBD&                a_dest,
                     const Copier&       a_copier,
                     const LDOP&         a_op = LDOP()) const; 

  void makeItSoLocalCopy(const LBD&          a_src,
                         LBD&                a_dest,
                         const Copier&       a_copier,
                         const LDOP&         a_op = LDOP()) const;

  void makeItSoEnd(const LDOP&         a_op = LDOP());
  
  //========================================================================
  //
  // data structures used by makeItSo when we have some
  // data that needs to be moved (ie. there are entries
  // in the 'FROM' or 'TO' CopyIterators)
  //
  void completePendingSends() const;

  void allocateBuffers(const LBD& a_src, 
                       const LBD& a_dest,
                       const Copier&   a_copier,
                       const LDOP& a_op) const;

  void writeSendDataFromMeIntoBuffers(const LBD& a_src,
                                      const LDOP& a_op) const;
  void postSendsFromMe() const ;

  void postReceivesToMe() const ;

  void unpackReceivesToMe(const LDOP& a_op) ;

  mutable shared_ptr<CopierBuffer> m_buff;
  };
#include "implem/Proto_LevelBoxDataImplem.H"
} //end namespace Proto
#endif
