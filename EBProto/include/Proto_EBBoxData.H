
#ifndef _Proto_EBBoxDATA_H_
#define _Proto_EBBoxDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{
  struct EBDataLoc
  {
    int    m_dataID;   //which pointer to index
    size_t m_offset;   //distance from start of pointer 
  };

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class BaseEBBoxData 
  {
  public:
    static const bool s_stackDefault = false;

    ///
    BaseEBBoxData()
    {
    }

    virtual ~BaseEBBoxData()
    {
    }

    ///
    BaseEBBoxData<cent, data_t, ncomp>&
    operator-=(const BaseEBBoxData<cent, data_t, ncomp>& a_input)
    {
      m_regData -= a_input.m_regData;
      m_irrData -= a_input.m_irrData;

      return *this;
    }


    ///
    BaseEBBoxData<cent, data_t, ncomp>&
    operator+=(const BaseEBBoxData<cent, data_t, ncomp>& a_input)
    {
      m_regData += a_input.m_regData;
      m_irrData += a_input.m_irrData;

      return *this;
    }
    
    ///
    BaseEBBoxData<cent, data_t, ncomp>&
    operator*=(const BaseEBBoxData<cent, data_t, ncomp>& a_input)
    {
      m_regData *= a_input.m_regData;
      m_irrData *= a_input.m_irrData;

      return *this;
    }


    ///
    virtual  void setVal(data_t a_val)
    {
      m_regData.setVal(a_val);
      m_irrData.setVal(a_val);
    }
    



    
    ///will not work on device--debugging use only
    inline  data_t operator()(const EBIndex<cent>& a_vof, 
                              unsigned int  a_c) const
    {
      if(m_graph.isIrregular(a_vof.m_pt))
      {
        return m_irrData(a_vof, a_c);
      }
      else
      {
        return m_regData(a_vof.m_pt, a_c);
      }
    }

    ///for AggStencil
    inline const data_t* dataPtr(int a_dataType, int a_ivar) const
    {
      const data_t*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
      return retval;
    }


    ///for AggStencil
    inline data_t* dataPtr(int a_dataType, int a_ivar) 
    {
      data_t*  retval = nullptr;
      if(a_dataType == 0)
      {
        retval = m_regData.dataPtr(a_ivar);
      }
      else
      {
        retval = m_irrData.dataPtr(a_ivar);
      }
      return retval;
    }

    ///needed to apply stencils
    BoxData<data_t,ncomp>& getRegData() 
    {
      return m_regData;
    }

    ///needed to apply stencils
    const BoxData<data_t,ncomp>& getRegData() const
    {
      return m_regData;
    }


    ///needed to apply stencils
    IrregData<cent, data_t,ncomp>& getIrregData() 
    {
      return m_irrData;
    }

    ///needed to apply stencils
    const IrregData<cent, data_t, ncomp>& getIrregData() const
    {
      return m_irrData;
    }
    ///
    const PointSet& irreg() const
    {
      return m_irreg;
    }

    ///needed for Proto::LevelData.  meaningless here becasue we need the graph.
    BaseEBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }


    //aliasing define
    template<unsigned int srccomp>
    void define(BaseEBBoxData<cent, data_t, srccomp>& a_input,
                unsigned int                          a_comp)
    {
      m_inputBox = a_input.getRegData().box();
      m_graph = a_input.ebgraph();
      m_irreg = a_input.irreg();
      m_irrData.define(a_input.getIrregData(), a_comp);
      m_regData.define(a_input.getRegData()  , a_comp);
    }           


    
    ///for data transfers.  datasize does NOT include sizeof(data_t) (to conform with boxdata::size)
    void dataInfo(data_t*& a_dataPtr, size_t& a_datasize, unsigned int a_whichfab)
    {
      PR_assert((a_whichfab == 0) || (a_whichfab == 1));
      if(a_whichfab == 0)
      {
        a_dataPtr  =  m_regData.data();
        a_datasize =  m_regData.size();
      }
      else 
      {
        a_dataPtr  =  m_irrData.data();
        a_datasize =  m_irrData.size();
      }
    }

    /// input box is cell centered---semantic is somewhat different depending on centering
//    inline virtual  void define(const Box& a_box, const EBGraph& a_graph) = 0;

    Box box() const
    {
      return m_regData.box();
    }

    EBGraph ebgraph() const
    {
      return m_graph;
    }

    Box inputBox() const
    {
      return m_inputBox;
    }
  protected:
    BaseEBBoxData(const            BaseEBBoxData<cent, data_t, ncomp>& a_src);
    BaseEBBoxData& operator=(const BaseEBBoxData<cent, data_t, ncomp>& a_src);

    Box m_inputBox;
    BoxData<        data_t, ncomp>     m_regData;
    IrregData<cent, data_t, ncomp>     m_irrData;
    EBGraph                            m_graph;
    PointSet                           m_irreg;
  };


  ///Holder for EB data over a box
  /**
     This is the general holder for all EB.
     if cent == CELL, this is over cells.
     if cent == BOUNDARY, the data only exists at cut cells
     if cent == XFACE, YFACE, ZFACE, the data is at coordinate face.
     All box inputs are assumed to be cll centered
   */
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  class EBBoxData: public BaseEBBoxData<cent, data_t, ncomp>
  {
  public:
    ///
    EBBoxData()
    {
    }

    virtual ~EBBoxData()
    {
    }

    ///
    /**
       dummy constructor used by leveldata
     */
     EBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }

    /// 
    /**
         Input box (a_box)  is cell centered. the box is adjusted for centering if it is cent==XFACE YFACE or ZFACE.
         All irregular data (data at cut cells is held in  a sparse holder (IrregData).    Regular data is held in a BoxData.
    */
    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {

      this->m_inputBox = a_box;
      define(a_box, a_graph, a_useStack);
    }
    
    ///needed for leveldata
    EBBoxData(const Box& a_box)
    {  }

  

    /// 
    /**
         Input box (a_box)  is cell centered. the box is adjusted for centering if it is cent==XFACE YFACE or ZFACE.
         All irregular data (data at cut cells is held in  a sparse holder (IrregData).    Regular data is held in a BoxData.
       
     */
    inline void define(const Box& a_box, const EBGraph& a_graph, bool a_useStack= false);

    /// for AggStencil -- not really part of the public interface.
    inline EBDataLoc indexInfo(const EBIndex<cent>& a_baseInd, int a_ivar) const;

      
  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<CELL, data_t, ncomp>: public BaseEBBoxData<CELL, data_t, ncomp>
  {
  public:


    //aliasing define
    template<unsigned int srccomp>
    void define(EBBoxData<CELL, data_t, srccomp>& a_input,
                unsigned int                      a_comp)
    {
      BaseEBBoxData<CELL, data_t, ncomp>::define(a_input, a_comp);
    }
    
    EBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }


    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_stackAlloc = true)
    {
      define(a_box, a_graph, a_stackAlloc);
    }

    void
    define(const Box& a_box, const EBGraph& a_graph, bool a_stackAlloc = true)
    {
      this->m_inputBox = a_box;
      this->m_graph = a_graph;
      this->m_irreg = a_graph.getIrregCells(a_box);  

      this->m_irrData.define(a_box, this->m_graph);  
      this->m_regData.define(a_box, a_stackAlloc);  //stack allocation turned off for now
    }


    EBBoxData(const Box& a_box)
    {  }


    EBBoxData()
    {  }

    EBDataLoc
    indexInfo(const EBIndex<CELL>& a_index, int a_ivar) const
    {
      //this ought to work even if the index is outside the domain
      //should just send it to reg fab
      EBDataLoc retval;
      if(this->m_irreg.contains(a_index.m_pt))
      {
        retval.m_dataID   = 1;
        retval.m_offset   = this->m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = this->m_regData.index(a_index.m_pt, a_ivar);
      }
      return retval;
    }


  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<BOUNDARY, data_t, ncomp>: public BaseEBBoxData<BOUNDARY, data_t, ncomp>
  {
  public:

    EBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }


    EBBoxData()
    {  }


    //aliasing define
    template<unsigned int srccomp>
    void define(EBBoxData<BOUNDARY, data_t, srccomp>& a_input,
                unsigned int                      a_comp)
    {
      BaseEBBoxData<BOUNDARY, data_t, ncomp>::define(a_input, a_comp);
    }

    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_stackAlloc=true)
    {
      define(a_box, a_graph, a_stackAlloc);
    }

    EBBoxData(const Box& a_box)
    {  }


    void
    define(const Box& a_box, const EBGraph& a_graph, bool a_stackAlloc = true)
    {
      this->m_inputBox = a_box;
      this->m_graph = a_graph;
      this->m_irreg = a_graph.getIrregCells(a_box);  
      this->m_irrData.define(a_box, this->m_graph, a_stackAlloc);  

    }

    EBDataLoc
    indexInfo(const EBIndex<BOUNDARY>& a_index, int a_ivar) const
    {
      EBDataLoc retval;
      retval.m_dataID   = 1;
      retval.m_offset   = this->m_irrData.index(a_index, a_ivar);
      return retval;
    }

  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<XFACE, data_t, ncomp>: public BaseEBBoxData<XFACE, data_t, ncomp>
  {
  public:

    EBBoxData()
    {  }


    //aliasing define
    template<unsigned int srccomp>
    void define(EBBoxData<XFACE, data_t, srccomp>& a_input,
                unsigned int                      a_comp)
    {
      BaseEBBoxData<XFACE, data_t, ncomp>::define(a_input, a_comp);
    }

    EBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }

    /// input box is cell centered
    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {
      define(a_box, a_graph, a_useStack);
    }


    EBBoxData(const Box& a_box)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {
      this->m_inputBox = a_box;
      Box facebox = a_box.growHi(0, 1);
      this->m_graph = a_graph;
      this->m_irreg = a_graph.getIrregCells(a_box);  

      this->m_irrData.define(a_box, this->m_graph);  
      this->m_regData.define(facebox, this->s_stackDefault);  //stack allocation turned off for now

    }

    EBDataLoc
    indexInfo(const EBIndex<XFACE>& a_index, int a_ivar) const
    {
      Point hipt = a_index.m_pt;
      Point lopt = a_index.m_pt - Point::Basis(0, 1);
      EBDataLoc retval;
      if(this->m_irreg.contains(hipt) || this->m_irreg.contains(lopt) )
      {
        retval.m_dataID   = 1;
        retval.m_offset   = this->m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = this->m_regData.index(hipt, a_ivar);
      }
      return retval;
    }

  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<YFACE, data_t, ncomp>: public BaseEBBoxData<YFACE, data_t, ncomp>
  {
  public:

    EBBoxData()
    {  }


    //aliasing define
    template<unsigned int srccomp>
    void define(EBBoxData<YFACE, data_t, srccomp>& a_input,
                unsigned int                       a_comp)
    {
      BaseEBBoxData<YFACE, data_t, ncomp>::define(a_input, a_comp);
    }

    EBBoxData(const Box& a_box, bool a_stackAlloc, bool a_useStack = false)
    {
    }


    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {
      define(a_box, a_graph, a_useStack);
    }


    EBBoxData(const Box& a_box, bool a_useStack = false)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph, bool a_useStack = false)
    {
      this->m_inputBox = a_box;
      Box facebox = a_box.growHi(1, 1);
      this->m_graph = a_graph;
      this->m_irreg = a_graph.getIrregCells(a_box);  

      this->m_irrData.define(a_box, this->m_graph);  
      this->m_regData.define(facebox, false);  //stack allocation turned off for now

    }

    EBDataLoc
    indexInfo(const EBIndex<YFACE>& a_index, int a_ivar) const
    {
      Point hipt = a_index.m_pt;
      Point lopt = a_index.m_pt - Point::Basis(1, 1);
      EBDataLoc retval;
      if(this->m_irreg.contains(hipt) || this->m_irreg.contains(lopt) )
      {
        retval.m_dataID   = 1;
        retval.m_offset   = this->m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = this->m_regData.index(hipt, a_ivar);
      }
      return retval;
    }

  };


  template <typename  data_t, unsigned int ncomp>
  class EBBoxData<ZFACE, data_t, ncomp>: public BaseEBBoxData<ZFACE, data_t, ncomp>
  {
  public:

    EBBoxData(const Box& a_box, bool a_stackAlloc)
    {
    }


    //aliasing define
    template<unsigned int srccomp>
    void define(EBBoxData<ZFACE, data_t, srccomp>& a_input,
                unsigned int                       a_comp)
    {
      BaseEBBoxData<ZFACE, data_t, ncomp>::define(a_input, a_comp);
    }

    EBBoxData(const Box& a_box, const EBGraph& a_graph, bool a_useStack=false)
    {
      define(a_box, a_graph, a_useStack);
    }


    EBBoxData(const Box& a_box)
    {  }

    void
    define(const Box& a_box, const EBGraph& a_graph, bool a_useStack= false)
    {
      
#if DIM==3
      this->m_inputBox = a_box;
      Box facebox = a_box.growHi(2, 1);
      this->m_graph = a_graph;
      this->m_irreg = a_graph.getIrregCells(a_box);  

      this->m_irrData.define(a_box, this->m_graph);  
      this->m_regData.define(facebox, false);  //stack allocation turned off for now
#endif
    }

    EBDataLoc
    indexInfo(const EBIndex<ZFACE>& a_index, int a_ivar) const
    {
      Point hipt = a_index.m_pt;
      Point lopt = a_index.m_pt - Point::Basis(2, 1);
      EBDataLoc retval;
      if(this->m_irreg.contains(hipt) || this->m_irreg.contains(lopt) )
      {
        retval.m_dataID = 1;
        retval.m_offset = this->m_irrData.index(a_index, a_ivar);
      }
      else
      {
        retval.m_dataID   = 0;
        retval.m_offset   = this->m_regData.index(hipt, a_ivar);
      }
      return retval;
    }
  };


  template<typename T>
  inline T&
  getIrregData(T& a_s)
  {
    return a_s;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline IrregData<cent, data_t, ncomp>&
  getIrregData(EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getIrregData();
  }



  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline IrregData<cent, data_t, ncomp>&
  getIrregData(const EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return (IrregData<cent, data_t, ncomp>&)(a_s.getIrregData());
  }


  template<CENTERING cent,  typename data_t, unsigned int ncomp>
  struct
  EBIrregStruct
  {
    data_t*        m_startPtr;
    unsigned int   m_varsize;
    unsigned int   m_offset;
    Point          m_index;
  };

  template <CENTERING cent, typename T>
  inline T
  getEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                T& a_T)
  {
    return a_T;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline vector< EBIrregStruct<cent, data_t, ncomp> >
  getEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                   IrregData<cent, data_t, ncomp>& a_s )
  {
//    data_t*      debPtr  = a_s.data();
//  printf("irreg data ptr = %p\n", debPtr);
    vector< EBIrregStruct<cent, data_t, ncomp> > retval;
    for(int ivec = 0; ivec < a_indices.size(); ivec++)
    {
////begin debug
//      Point pt = a_indices[ivec].m_pt;
//      if((pt[0]==17) && pt[1]==18)
//      {
//        std::cout << "pt = " << pt << ", ivec = " << ivec << std::endl;
//      }
////end debug
      EBIrregStruct<cent, data_t, ncomp>  vecval;
      vecval.m_startPtr = a_s.data();
      vecval.m_varsize  = a_s.vecsize();
      vecval.m_offset   = a_s.index(a_indices[ivec], 0);
      vecval.m_index    = a_indices[ivec].m_pt;
      retval.push_back(vecval);
    }

//  printf("host return data ptr = %p\n", retval.data());
    return retval;
  }


  template<typename T>
  inline T&
  getBoxData(T& a_s)
  {
    return a_s;
  }

  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline BoxData<data_t, ncomp>&
  getBoxData(EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getRegData();
  }


  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline const BoxData<data_t, ncomp>&
  getBoxData(const EBBoxData<cent, data_t, ncomp>& a_s)
  {
    return a_s.getRegData();
  }


#ifdef PROTO_CUDA

  template<typename T>
  __device__ __host__
  inline T
  cudaGetVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }


  template<CENTERING cent, typename data_t, unsigned int ncomp>
  __device__ 
  inline Var<data_t, ncomp>
  cudaGetVar(unsigned int a_ivec,
             EBIrregStruct<cent, data_t, ncomp>* a_dst)
  {
    Var<data_t, ncomp> retval;

    const EBIrregStruct<cent, data_t, ncomp>*  rawptr = a_dst;
    const EBIrregStruct<cent, data_t, ncomp>&  ugly   = rawptr[a_ivec];
//  printf("cudaGetVar rawptr   = %p \n", rawptr);
//  printf("cudaGetVar startptr = %p \n", ugly.m_startPtr);
    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
    }
    return retval;
  }




//going into this srcs are thrust_device_pointer<uglystruct> and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  __global__
  void
  vec_indexer(unsigned int a_begin, unsigned int a_end,Func a_body, 
              EBIrregStruct<cent, data_t, ncomp>*  a_dst, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
    {
      a_body(cudaGetVar(idx, a_dst), cudaGetVar(idx, a_srcs)...);
    }
  }


  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  __global__
  void
  vec_indexer_i(unsigned int a_begin, unsigned int a_end,Func a_body, 
                EBIrregStruct<cent, data_t, ncomp>*  a_dst, Srcs... a_srcs)
  {
    int idx = threadIdx.x + blockIdx.x*blockDim.x;
    if (idx >= a_begin && idx < a_end)
    {
      a_body(a_dst[idx].m_index.m_tuple, cudaGetVar(idx, a_dst), cudaGetVar(idx, a_srcs)...);
    }
  }
//going into this srcs are uglystruct* and other stuff
  template <typename T>
  inline int
  cleanUpPtrMem(T a_T)
  {
    return 0;
  }
//going into this srcs are uglystruct* and other stuff
  template <CENTERING cent, typename data_t,unsigned int ncomp>
  inline int
  cleanUpPtrMem(EBIrregStruct<cent, data_t, ncomp> * a_ptr)
  {
    cudaFree(a_ptr);
    return 0;
  }

//
  template<typename... Srcs>
  inline void
  emptyFunc(Srcs... a_srcs)
  {
  }

//going into this srcs are uglystruct* and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  void
  cudaVectorFunc(const Func& a_F, unsigned int a_Nvec, 
                 EBIrregStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
  {
//  printf("cudavecf: dst  = %p\n", a_dst);
    //printf("cudavecf: src  = %p\n", a_firstsrc);
    cudaStream_t curstream = DisjointBoxLayout::getCurrentStream();
    const int N = a_Nvec;
    unsigned int stride = a_Nvec;
    unsigned int blocks = 1;
    size_t smem = 0;
    vec_indexer<<<blocks, stride, smem, curstream>>>
      (0, N, mapper(a_F), a_dst, a_srcs...);

    //there is a cudaMalloc that happens above so we have to delete
    emptyFunc(cleanUpPtrMem(a_dst ), (cleanUpPtrMem(a_srcs))...); 
  }


//going into this srcs are uglystruct* and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  void
  cudaVectorFunc_i(const Func& a_F, unsigned int a_Nvec, 
                   EBIrregStruct<cent, data_t, ncomp> * a_dst,Srcs... a_srcs)
  {
//  printf("cudavecf: dst  = %p\n", a_dst);
    //printf("cudavecf: src  = %p\n", a_firstsrc);
    cudaStream_t curstream = DisjointBoxLayout::getCurrentStream();
    const int N = a_Nvec;
    unsigned int stride = a_Nvec;
    unsigned int blocks = 1;
    size_t smem = 0;
    vec_indexer_i<<<blocks, stride, smem, curstream>>>
      (0, N, mapper(a_F), a_dst, a_srcs...);

    //there is a cudaMalloc that happens above so we have to delete
    emptyFunc(cleanUpPtrMem(a_dst ), (cleanUpPtrMem(a_srcs))...); 
  }


  template <CENTERING cent, typename T>
  inline T
  cudagetEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                    T& a_T)
  {
    return a_T;
  }
//
  template <CENTERING cent, typename  data_t, unsigned int ncomp>
  inline  EBIrregStruct<cent, data_t, ncomp>*
  cudagetEBIrregStruct(const vector<EBIndex<cent> >& a_indices,
                    IrregData<cent, data_t, ncomp>& a_s )
  {
    vector< EBIrregStruct<cent, data_t, ncomp> > hostvec = getEBIrregStruct(a_indices, a_s);

    size_t memsize = hostvec.size()*sizeof(EBIrregStruct<cent, data_t, ncomp>);
    EBIrregStruct<cent, data_t, ncomp>* retval;
    cudaMalloc(&retval, memsize);
    cudaMemcpy(retval, hostvec.data(), memsize, cudaMemcpyHostToDevice);

    //this copies from the host to the device
//  printf("cgus: device host vector ptr = %p\n", hostvec.data());
//  printf("cgus: device return data ptr = %p\n", retval);
    return retval;
  }
//going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  cudaEBForAllIrreg(const Func& a_F, const Box& a_box,
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
  {
    //indicies into irreg vector that correspond to input box
    const vector<EBIndex<cent> >& dstvofs = a_dst.getIndices(a_box);
    unsigned int vecsize = a_dst.vecsize();
    if(vecsize > 0)
    {
//    printf("cudaebforall: dst  = %p\n", a_dst.data());
      cudaVectorFunc(a_F, vecsize, cudagetEBIrregStruct(dstvofs, a_dst), 
                     (cudagetEBIrregStruct(dstvofs, a_srcs))...);

    }
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  cudaEBForAllIrreg_i(const Func& a_F, const Box& a_box,
                      IrregData<cent, data_t, ncomp>& a_dst,
                      Srcs&...  a_srcs)
  {
    //indicies into irreg vector that correspond to input box
    const vector<EBIndex<cent> >& dstvofs = a_dst.getIndices(a_box);
    unsigned int vecsize = a_dst.vecsize();
    if(vecsize > 0)
    {
//    printf("cudaebforall: dst  = %p\n", a_dst.data());
      cudaVectorFunc_i(a_F, vecsize, cudagetEBIrregStruct(dstvofs, a_dst), 
                       (cudagetEBIrregStruct(dstvofs, a_srcs))...);

    }
  }

  ///
  template<typename Func, typename... Srcs>
  inline void
  cudaEBforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    cudaEBForAllIrreg(a_F, a_box, getIrregData(a_srcs)...);
  }


  template<typename Func, typename... Srcs>
  inline void
  cudaEBforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    cudaEBForAllIrreg_i(a_F, a_box, getIrregData(a_srcs)...);
  }

#else
///cpu-only specific functions

///
  template<typename T>
  inline T
  getVar(unsigned int ivec,  T a_s)
  {
    return a_s;
  }

///
  template<CENTERING cent, typename data_t, unsigned int ncomp>
  inline Var<data_t, ncomp>
  getVar(unsigned int a_ivec,
         vector< EBIrregStruct<cent, data_t, ncomp> > a_dst)
  {
    Var<data_t, ncomp> retval;
    const EBIrregStruct<cent, data_t, ncomp> ugly = a_dst[a_ivec];
    for(int icomp = 0; icomp < ncomp; icomp++)
    {
      retval.m_ptrs[icomp] = ugly.m_startPtr + ugly.m_offset + (ugly.m_varsize*icomp);
    }
//begin debug
//    if(a_ivec == 19)
//    {
//      std::cout << "getvar, ivec == 19, values = "  ;
//      for(int icomp = 0; icomp < ncomp; icomp++)
//      {
//        std::cout << *retval.m_ptrs[icomp] << " ";
//      }
//      std::cout << std::endl;
//    }
//end debug

    return retval;
  }

///going into this srcs are vector<uglystruct> and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  void
  hostVectorFunc(const Func& a_F, vector< EBIrregStruct<cent, data_t, ncomp> > a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
    {
      a_F(getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }


///going into this srcs are vector<uglystruct> and other stuff
  template<CENTERING cent, typename data_t,unsigned int ncomp,  typename Func, typename... Srcs>
  void
  hostVectorFunc_i(const Func& a_F, vector< EBIrregStruct<cent, data_t, ncomp> > a_dst, Srcs... a_srcs)
  {
    for(unsigned int ivec = 0; ivec < a_dst.size(); ivec++)
    {
      Point pt = a_dst[ivec].m_index;
      a_F(pt.m_tuple, getVar(ivec, a_dst), (getVar(ivec, a_srcs))...);
    }
       
  }

///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBForAllIrreg(const Func& a_F, const Box& a_box,
                    IrregData<cent, data_t, ncomp>& a_dst,
                    Srcs&...  a_srcs)
  {
//indicies into irreg vector that correspond to input box
    vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
    hostVectorFunc(a_F, getEBIrregStruct(dstvofs, a_dst), (getEBIrregStruct(dstvofs, a_srcs))...);
  }


///going into this srcs are IrregDatas and other stuff
  template<CENTERING cent, typename  data_t, unsigned int ncomp, typename Func, typename... Srcs>
  inline void
  hostEBForAllIrreg_i(const Func& a_F, const Box& a_box,
                      IrregData<cent, data_t, ncomp>& a_dst,
                      Srcs&...  a_srcs)
  {
//indicies into irreg vector that correspond to input box
    vector<EBIndex<cent> > dstvofs = a_dst.getIndices(a_box);
    hostVectorFunc_i(a_F, getEBIrregStruct(dstvofs, a_dst), (getEBIrregStruct(dstvofs, a_srcs))...);
  }


///going into this srcs are EBBoxDatas and other stuff
  template<typename Func, typename... Srcs>
  inline void
  hostEBforall(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
    forallInPlaceBase(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
    hostEBForAllIrreg(a_F, a_box, (getIrregData(a_srcs))...);
  }

///going into this srcs are EBBoxDatas and other stuff
  template<typename Func, typename... Srcs>
  inline void
  hostEBforall_i(const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
//call regular forall
//    std::cout << "doing regular   forall_i" << std::endl;
    forallInPlace_i(a_F, a_box, (getBoxData(a_srcs))...);
  
//do the same thing for the irregular data
//    std::cout << "doing irregular forall_i" << std::endl;
    hostEBForAllIrreg_i(a_F, a_box, getIrregData(a_srcs)...);
  }

#endif

///version that does not send the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace(unsigned long long int a_num_flops_point,
                              const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

//  printf("in ebforall function pointer = %p\n", &a_F);
    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
    cudaEBforall(a_F, a_box, a_srcs...);
    cudaDeviceSynchronize();
#else
    hostEBforall(a_F, a_box, a_srcs...);
#endif
    PR_FLOPS(boxfloops);
  }


///version that sends the point to the function
  template<typename Func, typename... Srcs>
  inline void ebforallInPlace_i(unsigned long long int a_num_flops_point,
                                const char*            a_timername,
                                const Func & a_F,  Box a_box, Srcs&... a_srcs)
  {
    PR_TIME(a_timername);

    unsigned long long int boxfloops = a_num_flops_point*a_box.size();

#ifdef PROTO_CUDA
    cudaEBforall_i(a_F, a_box, a_srcs...);
#else
    hostEBforall_i(a_F, a_box, a_srcs...);
#endif


    PR_FLOPS(boxfloops);
  }

  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg(const char*            a_timername,
                            const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
#ifdef PROTO_CUDA
    cudaEBForAllIrreg(a_F, a_box, getIrregData(a_srcs)...);
#else
    hostEBForAllIrreg(a_F, a_box, (getIrregData(a_srcs))...);
#endif
  }



  ///version that only touches cut cells
  template<typename Func, typename... Srcs>
  inline void ebforallIrreg_i(const char*            a_timername,
                              const Func & a_F,  Box a_box, Srcs&... a_srcs)

  {
    PR_TIME(a_timername);
    //the flops thing here is tricky
#ifdef PROTO_CUDA
    cudaEBForAllIrreg_i(a_F, a_box, getIrregData(a_srcs)...);
#else
    hostEBForAllIrreg_i(a_F, a_box, (getIrregData(a_srcs))...);
#endif
  }
}

#endif

