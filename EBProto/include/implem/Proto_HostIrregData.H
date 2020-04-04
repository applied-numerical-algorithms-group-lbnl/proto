
#ifndef _HostIrregData_H_
#define _HostIrregData_H_

using std::vector;
using std::shared_ptr;
#include "EBProto.H"
#include "Proto_EBGraph.H"
#include "Proto_Serialization.H"
namespace Proto
{
///
/**
*/
///
/**
 */
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class HostIrregData
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    HostIrregData()
    {
      m_boxdefined = false;
    }

    ///
    virtual ~HostIrregData()
    {
    }

    ///
    HostIrregData( const shared_ptr< vector<index_t> >& a_indicies)
    {
      define(a_indicies);
    }

    ///needed for leveldata
    HostIrregData(const Box& a_box, bool a_stackAllocation=true)
    {  }

    ///
    inline void define( const shared_ptr< vector< index_t > >& a_indicies)
    {
      m_indicies = a_indicies;
      std::size_t count = this->size();
      m_data=::std::shared_ptr<data_t>(new data_t[count], [](data_t* p) { delete[] p;});
      m_rawPtr = m_data.get();
    }


    ///
    inline std::size_t size() const {return m_indicies->size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    data_t* dataPtr(unsigned int a_ivar = 0)
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    ///
    const data_t* dataPtr(unsigned int a_ivar = 0) const
    {
      PR_assert(a_ivar < ncomp);
      size_t ioff = a_ivar*m_indicies->size();
      return &(m_rawPtr[ioff]);
    }

    data_t* data()
    {
      return m_rawPtr;
    }


    const data_t* data() const
    {
      return m_rawPtr;
    }
    
    ///
    inline size_t  index(const EBIndex<cent>& a_vof, const int& a_ivar) const
    {
      size_t ivec = 0;
      bool found = false;
      for(unsigned int i = 0; i < m_indicies->size(); i++)
      {        
        if(a_vof == (*m_indicies)[i])
        {
          ivec = i;
          found = true;
          break;
        }
      }
      if(!found)
      {
        std::cout << "index<cell> vof not found" << std::endl;
        abort();
      }
      PROTO_ASSERT(found, "vof not found");

      ivec += a_ivar*m_indicies->size();

      if(ivec >= size())
      {
        std::cout << "index<cell> offset too big" << std::endl;
        abort();
      }

      return ivec;
    }

    /// 
    const data_t& operator() (const index_t& a_vof, int a_comp) const
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }
    
    data_t& operator() (const index_t& a_vof, int a_comp)
    {
      size_t ind = index(a_vof, a_comp);
      return m_rawPtr[ind];
    }

    ///might need a different semantic for face-centered data.
    vector<EBIndex<cent> > getIndices(const Box& a_box) const
    {
      vector<EBIndex<cent> > retval;
      for(int i = 0; i < m_indicies->size(); i++)
      {
        if(a_box.contains((*m_indicies)[i].m_pt))
        {
          retval.push_back((*m_indicies)[i]);
        }
      }
      return retval;
    }

    unsigned int vecsize() const
    {
      return m_indicies->size();
    }

    //aliasing define
    template<unsigned int srccomp>
    void define(HostIrregData<cent, data_t, srccomp>& a_input,
                unsigned int                        & a_comp)
    {
      m_indicies = a_input.getIndicies();
      m_data     = a_input.getData();
      m_rawPtr   = a_input.dataPtr(a_comp);
    }           

    //do not use this
    shared_ptr<vector<  index_t > >  getIndicies() const
    {
      return m_indicies;
    }

    //do not use this
    shared_ptr<data_t> getData() const
    {
      return  m_data;
    }

    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      this->m_boxdefined = true;
      this->m_inputBox   = a_box;
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector<EBIndex<cent> > > indicies(new vector< EBIndex<cent> >() );
      *indicies = getIrregIndicies<cent>(interBox, a_graph);

      define(indicies);
    }

    inline bool hasIndex(const EBIndex<cent>& a_ind) const
    {
      bool found = false;
      //this should probably be a smarter search.   maybe sort them?
      for(int ivec = 0;(ivec < m_indicies->size() && !found); ivec++)
      {
        if((*m_indicies)[ivec] == a_ind)
        {
          found = true;
        }
      }

      return found;
    }
    
    bool contains(const EBIndex<cent>& a_index) const
    {
      return hasIndex(a_index);
    }

    void setVal(data_t val)
    {
      for(int ivec = 0;ivec < m_indicies->size() ; ivec++)
      {
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
          (*this)((*m_indicies)[ivec], icomp) = val;
        }
      }
    }

    bool copyContains(const index_t& a_ind, const Box & a_srcbox) const
    {
      CENTERING thiscent = cent;

      bool retval = false;
      if((thiscent == CELL) || (thiscent == BOUNDARY))
      {        
        if(a_srcbox.contains(a_ind.m_pt))
        {
          retval = this->contains(a_ind);
        }
      }
      else if(thiscent == XFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(0, 1);
        if(a_srcbox.contains(hipt) || a_srcbox.contains(lopt))
        {
          retval = this->contains(a_ind);
        }
      }
      else if(thiscent ==  YFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(1, 1);
        if(a_srcbox.contains(lopt) || a_srcbox.contains(hipt))
        {
          retval = this->contains(a_ind);
        }
      }
#if DIM==3      
      else if(thiscent ==  ZFACE)
      {
        Point hipt = a_ind.m_pt;
        Point lopt = a_ind.m_pt - Point::Basis(2, 1);
        if(a_srcbox.contains(lopt) || a_srcbox.contains(hipt))
        {
          retval = this->contains(a_ind);
        }
      }
#endif
      else
      {
        PROTO_ASSERT(false, "unaccounted-for centering");
      }
      return retval;
    }

    ///this DOES include sizeof(T)  
      std::size_t charsize(const Box&   a_bx,
                           unsigned int a_startcomp,
                           unsigned int a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      using ProtoSer::linearCharBufSize;
      unsigned int numpts = 0;
      size_t retval = sizeof(unsigned int); //for the number of points transferred
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        const index_t & ind = (*m_indicies)[i];
        if(this->copyContains(ind, a_bx))
          {
            numpts++;
            for(int icomp = 0; icomp < a_numcomps; icomp++)
            {
              const data_t& dataval = (*this)(ind, a_startcomp + icomp);
              retval += linearCharBufSize<data_t>(dataval);
            }
          } 
      }
      retval += numpts*sizeof(index_t); // for the indices themselves

//      std::cout << "hostirregdata::charsize bx ="  << a_bx << ". size = " << retval << std::endl;

      return retval;
    }

    void linearOut(void*               a_buf,
                   const Box       &   a_bx,
                   unsigned int        a_startcomp,
                   unsigned int        a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      using ProtoSer::linearCharBufSize;
      using ProtoSer::linearOutToCharBuf;

      vector<index_t> indices;

      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        const index_t & ind = (*m_indicies)[i];
        if(this->copyContains(ind, a_bx))
          {
            indices.push_back(ind);
          } 
      }
      //output the size
      char* charbuf = (char*) a_buf;
      unsigned int* intbuf = (unsigned int*) charbuf;
      *intbuf = indices.size();
      charbuf += sizeof(unsigned int);

      for(int ipt = 0; ipt < indices.size(); ipt++)
      {
        //output the index
        const index_t& index = indices[ipt];

        index_t* indptr = (index_t *) charbuf;
        *indptr = index;
        charbuf += sizeof(index_t);

        //output the data
        for(int icomp = 0; icomp < a_numcomps; icomp++)
        {
          const data_t& dataval = (*this)(index, a_startcomp + icomp);

          linearOutToCharBuf<data_t>(charbuf, dataval);
          charbuf += linearCharBufSize<data_t>(dataval);
        }
      }
    }

    void  linearIn(void*                  a_buf,
                   const Box       &      a_bx,
                   unsigned int           a_startcomp,
                   unsigned int           a_numcomps)
    {
      PR_assert((a_startcomp + a_numcomps-1) < ncomp);
      using ProtoSer::linearCharBufSize;
      using ProtoSer::linearInFromCharBuf;
      char* charbuf = (char *) a_buf;
      int* intbuf = (int*)charbuf;
      unsigned int numindex = *intbuf;
      charbuf += sizeof(unsigned int);//for the number of points transferred
      for(int ipt = 0; ipt < numindex; ipt++)
      {
        index_t* indptr = (index_t*) charbuf;
        index_t  index = *indptr;
        charbuf += sizeof(index_t);

        //have to take them out of the buffer in any case since they were put in.
        //after that we check to see if we have a place to put them
        std::vector<data_t> datavals(a_numcomps);
        for(int icomp = 0; icomp < a_numcomps; icomp++)
        {
          data_t* dataptr = (data_t*) charbuf;
          linearInFromCharBuf<data_t>(*dataptr, charbuf);
          charbuf += linearCharBufSize<data_t>(*dataptr);
        }
        if(this->contains(index))
        {
          for(int icomp = 0; icomp < a_numcomps; icomp++)
          {
            data_t& dataval = (*this)(index, a_startcomp + icomp);
            char* tempcharbuf = (char*)(&datavals[icomp]);
            linearInFromCharBuf<data_t>(dataval, tempcharbuf);  //might not have operator=
          }
        }
      }
    }
    ///
    static int preAllocatable()
    {
      return 1; 
    }
    
    inline
    void copy(const HostIrregData<cent, data_t, ncomp>& a_src,
              const Box& a_bxFrom, unsigned int a_isrc,
              const Box& a_bxTo  , unsigned int a_idst,
              unsigned int a_inco)
    {
      PR_assert("HostIrregData::copy");
      Box interBox = a_bxFrom & a_bxTo;

      for(int iind = 0; iind < m_indicies->size(); iind++)
      {
        auto& vof = (*m_indicies)[iind];
        if(a_src.contains(vof) && (interBox.contains(vof.m_pt)))
        {
          for(unsigned int ico = 0; ico < a_inco; ico++)
          {
            unsigned int isrc = a_isrc + ico;
            unsigned int idst = a_idst + ico;
            (*this)(vof, idst) = a_src(vof, isrc);
          }
        }
      }
    }
  protected:

    bool m_boxdefined;
    Box m_inputBox; //exists if defined as part of a boxdata
    shared_ptr<vector<  index_t > >        m_indicies;
    shared_ptr<data_t>         m_data;
    data_t*                    m_rawPtr;
    HostIrregData(const HostIrregData<cent, data_t,ncomp>& a_src);
    HostIrregData& operator=(const HostIrregData  <cent, data_t,ncomp>& a_src);


  };

}

#endif
