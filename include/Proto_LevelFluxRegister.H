#include "Proto.H"
#include "Proto_Register.H"
#include "Proto_BoxData.H"
#include "Proto_Copier.H"
#define LFR LevelFluxRegister< T , C , MEMTYPE>
#define LFRT template<typename T, unsigned int C, MemType MEMTYPE>
#define LFROP LFRCopierOp<T,C,MEMTYPE>
namespace Proto
{
  LFRT
  class LevelFluxRegister;
  LFRT
  class LFRCopierOp
  {
  public:
    LFR* m_lfr; 
    ~LFRCopierOp()
    {
    }
    inline LFRCopierOp(LFR& a_lfr)
    {
      m_lfr = &a_lfr;
    }
    inline int size(const Box& a_bx) const
    {
      return C*a_bx.size()*sizeof(T);
    }
    inline void linearOut(void* a_buf, 
                          const Box& a_bx,
                          const DataIndex& a_diFr) const
    {
      // linearOut copies from source data holders into buffers.
      BoxData<T,C,MEMTYPE>& bd = m_lfr->sourceDB(a_bx,a_diFr);
      CInterval cint(0,C-1);
      bd.linearOut(a_buf, a_bx,cint);
    }
    inline void linearIn(void* a_buf, 
                         const Box& a_bx,
                         const DataIndex& a_diTo) const
    {
      // linearIn copies from buffers into destination data holders.
      BoxData<T,C>& bd = m_lfr->destDB(a_bx,a_diTo);
      CInterval cint(0,C-1);
      bd.linearIn(a_buf, a_bx,cint);
    }

    inline void op(const Box& a_bxFr,const DataIndex& a_diFr,
                   const Box& a_bxTo,const DataIndex& a_diTo) const
    {
      // Copies from source to data directly assuming they are on the same processor.
      int LinearizationTest = 1;
      if (LinearizationTest == 0)
        {
          PR_TIMERS("LOP_localCopy_1");
          BoxData<T,C>& src = m_lfr->sourceDB(a_bxFr,a_diFr);
          BoxData<T,C>& dest =m_lfr->destDB(a_bxTo,a_diTo);
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
  template
  <typename T, unsigned int C=1, MemType MEMTYPE=HOST>
  class LevelFluxRegister
  {
  public:
    /// Constructor.
    LevelFluxRegister(const DisjointBoxLayout& a_dblCoarse,
                      const DisjointBoxLayout& a_dblFine,
                      const Point& a_refRatio);

    /// Define function.
    void define(const DisjointBoxLayout& a_dblCoarse,
                const DisjointBoxLayout& a_dblFine,
                const Point& a_refRatio);


    /// Increments all the registers contained in the coarse patch a_di, direction a_dir. 
    /** Assumes that a_flux is defined on all the a_dir faces in the patch.
     */
    void incrementCoarse(const BoxData<T,C,MEMTYPE>& a_flux,
                         const DataIndex& a_di,
                         const T& a_weight,
                         unsigned int a_dir);

    /// Increments all the registers contained in the fine patch a_di, direction a_dir. 
    /// Assumes that a_flux is defined on all the a_dir faces in the patch.
    void incrementFine(const BoxData<T,C,MEMTYPE>& a_flux,
                       const DataIndex& a_di,
                       const T& a_weight,
                       unsigned int a_dir); 

    /// Increments coarse data using flux registers.
    void reflux(LevelBoxData<T,C,MEMTYPE>& a_coarseData,
                const T& a_weight);
    /// access to fine register data. Used in LFROperator for Copier.
    BoxData<T,C,MEMTYPE>& sourceBD(const Box& a_bx,const DataIndex& di);
    /// Access to coarse register data. Used in LFROperator for Copier.
    BoxData<T,C,MEMTYPE>& destBD(const Box& a_bx,const DataIndex& di);
    void reset();
    /// compute the hash key for a coarse register.
    int key(const Box& a_bx,const DataIndex& a_di);
  protected:
    DisjointBoxLayout m_dblFine,m_dblCoarse;
    Copier m_lfrCopier;
    Point m_refRatio;
    vector<vector<Register<T,C,MEMTYPE> > > m_coarseRegisters;
    vector<vector<shared_ptr<BoxData<T,C,MEMTYPE> > > > m_destData;
    vector<vector<Register<T,C,MEMTYPE> > > m_fineRegisters; 
    vector<unordered_map<unsigned int,unsigned int> > m_coarseIndices;
    array<array<Stencil<T>,2>,DIM> m_avgFlux;
    vector<array<array<int,2>, DIM> > m_hasFineRegister;

 // Internal fucntions that execute copying, exchange.
 void makeItSo(
               const Copier&       a_copier,
               const LFROP&         a_op) const;
  
  void makeItSoBegin(
                     const Copier&       a_copier,
                     const LFROP&         a_op) const; 

  void makeItSoLocalCopy(
                         const Copier&       a_copier,
                         const LFROP&         a_op) const;

  void makeItSoEnd(const LFROP&         a_op) const;
  
  //========================================================================
  //
  // data structures used by makeItSo when we have some
  // data that needs to be moved (ie. there are entries
  // in the 'FROM' or 'TO' CopyIterators)
  //
  void completePendingSends() const;

  void allocateBuffers(
                       const Copier&   a_copier,
                       const LFROP& a_op) const;

  void writeSendDataFromMeIntoBuffers(const LFROP& a_op) const;
  void postSendsFromMe() const ;

  void postReceivesToMe() const ;

  void unpackReceivesToMe(const LFROP& a_op) const;

  mutable CopierBuffer* m_buff;
  };
#include "implem/Proto_LevelFluxRegisterImplem.H"
}
