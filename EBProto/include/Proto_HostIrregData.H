
#ifndef _HostIrregData_H_
#define _HostIrregData_H_

using std::vector;
using std::shared_ptr;
#include "EBProto.H"
#include "Proto_EBGraph.H"
namespace Proto
{
///
/**
*/
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class BaseHostIrregData 
  {
  public:
    typedef EBIndex<cent> index_t;

    ///
    virtual ~BaseHostIrregData()
    {
    }

    ///
    BaseHostIrregData()
    {
    }

    ///
    BaseHostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }


    ///needed for leveldata
    BaseHostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {

    }

 
    void copy(const BaseHostIrregData<cent, data_t, ncomp>   &  a_src,
              const Box                                      &  a_srcbox,
              unsigned int                                      a_srccomp,
              const Box                                      &  a_dstbox,
              unsigned int                                      a_dstcomp,
              unsigned int                                      a_numcomp)
    {
      PR_assert((a_srccomp + a_numcomp) < ncomp)
      PR_assert((a_dstcomp + a_numcomp) < ncomp)
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        const index_t & ind = (*m_indicies)[i];
        if((a_src.copyContains(ind, a_srcbox)) && (this->copyContains(ind, a_dstbox)))
          {
            for(unsigned int icomp = 0; icomp < a_numcomp; icomp++)
            {
              unsigned int isrc = a_srccomp + icomp;
              unsigned int idst = a_dstcomp + icomp;
              (*this)(ind, idst) = a_src(ind, isrc);
            }
          }
      }
    }

    //what to do with box is specific to the implementation
    virtual bool copyContains(const index_t& a_ind, const Box & a_srcbox) const = 0;

    virtual bool contains(const index_t& a_ind) const
    {
      bool retval = false;
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        if((*m_indicies)[i] == a_ind)
        {
          retval = true;
          break;
        }
      }
      return retval;
    }

    ///
    inline void define(const shared_ptr<vector<index_t> >& a_indicies)
    {
      m_indicies = a_indicies;
      std::size_t count = this->size();

      m_data=::std::shared_ptr<data_t>(new data_t [count], [](data_t* p) { delete[] p;});
      m_rawPtr = m_data.get();
    }

    ///
    inline std::size_t size() const 
    {
      return m_indicies->size()*ncomp;
    }

    ///
    shared_ptr< vector<index_t> > getIndicies() const
    {
      return m_indicies;
    }


    ///this DOES include sizeof(T)  
    virtual std::size_t charsize(const Box&   a_bx,
                                 unsigned int a_startcomp,
                                 unsigned int a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
      unsigned int numpts = 0;
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        const index_t & ind = (*m_indicies)[i];
        if(this->copyContains(ind, a_bx))
          {
            numpts++;
          } 
      }
      size_t retval = sizeof(unsigned int); //for the number of points transferred
      retval += numpts*sizeof(index_t); // for the indices themselves
      retval += numpts*sizeof(data_t)*a_numcomps; //for the actual data
      return retval;
    }

    virtual void linearOut(void*               a_buf,
                           const Box       &   a_bx,
                           unsigned int        a_startcomp,
                           unsigned int        a_numcomps) const
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
      vector<index_t> indices;

      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        const index_t & ind = (*m_indicies)[i];
        if(this->copyContains(ind, a_bx))
          {
            indices.push_back(ind);
          } 
      }
      char* charbuf = (char*) a_buf;
      unsigned int* intbuf = (unsigned int*) charbuf;
      *intbuf = indices.size();
      charbuf += sizeof(unsigned int);

      for(int ipt = 0; ipt < indices.size(); ipt++)
      {
        index_t* indptr = (index_t *) charbuf;
        index_t index = *indptr;
        charbuf += sizeof(index_t);
        for(int icomp = 0; icomp < a_numcomps; icomp++)
        {
          data_t* dataptr = (data_t*) charbuf;
          *dataptr = (*this)(index, a_startcomp + icomp);
          charbuf += sizeof(data_t);
        }
      }
    }

    virtual void  linearIn(void*                  a_buf,
                           const Box       &      a_bx,
                           unsigned int           a_startcomp,
                           unsigned int           a_numcomps)
    {
      PR_assert((a_startcomp + a_numcomps) < ncomp);
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
          datavals[icomp] = *dataptr;
          charbuf += sizeof(data_t);
        }
        if(this->contains(index))
        {
          for(int icomp = 0; icomp < a_numcomps; icomp++)
          {
            (*this)(index, a_startcomp + icomp) = datavals[icomp];
          }
        }
      }
    }


    ///
    data_t* data()
    {
      return m_rawPtr;
    }

    ///
    const data_t* data() const
    {
      return m_rawPtr;
    }

    /**
     */
    virtual data_t* getIndex(const index_t& a_vof,const int& a_comp) const
    {
      PR_assert(a_comp < ncomp);
      data_t* dataPtr =  NULL;
      bool found      = false;
      data_t* ptrloc   = m_rawPtr;
      for (unsigned int i=0; i<m_indicies->size(); ++i)
      {
        if (a_vof == (*m_indicies)[i]) 
        {
          found = true;
          dataPtr = ptrloc;
          break;
        }
        ptrloc++;
      }

      PROTO_ASSERT(found, "index not found for hostirregdata");

      if(found)
      {
        dataPtr = ptrloc + a_comp*m_indicies->size();
      }
      return dataPtr;
    }


    ///
    data_t& operator() (const index_t& a_vof, int a_comp)
    {  return *(getIndex(a_vof, a_comp));}

    ///
    const data_t& operator() (const index_t& a_vof, int a_comp) const
    {
      return *(getIndex(a_vof, a_comp));
    }

    ///
    inline void setVal(const data_t& a_val)
    {
      for(int ival = 0; ival < size(); ival++)
      {
        m_rawPtr[ival] = a_val;
      }
    }

  protected:

    BaseHostIrregData(const            BaseHostIrregData<cent, data_t, ncomp>& a_src);
    BaseHostIrregData& operator=(const BaseHostIrregData<cent, data_t, ncomp>& a_src);

    shared_ptr<vector<index_t> >  m_indicies;
    shared_ptr<data_t>            m_data;
    data_t*                       m_rawPtr;
  };


  ///
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class HostIrregData: public BaseHostIrregData<cent, data_t, ncomp> 
  {
  public:

    typedef EBIndex<cent> index_t;
    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    virtual bool copyContains(const index_t& a_ind, const Box & a_srcbox) const;

    HostIrregData()
    {
    }

    HostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph);
  };



  ///
  template <class data_t,  unsigned int ncomp>
  class HostIrregData<CELL, data_t, ncomp>: public BaseHostIrregData<CELL, data_t, ncomp> 
  {
  public:

    typedef EBIndex<CELL> index_t;


    virtual bool copyContains(const index_t& a_ind, const Box & a_box) const
    {

      bool retval = false;
      if(a_box.contains(a_ind.m_pt))
        {
          retval = this->contains(a_ind);
        }
      return retval;
    }

    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    HostIrregData()
    {
    }

    HostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector< EBIndex<CELL> > > indicies(new vector< EBIndex<CELL> >() );
      *indicies = a_graph.getIrregVoFs(interBox);

      BaseHostIrregData<CELL, data_t, ncomp>::define(indicies);
    }
  };


  ///
  template <class data_t,  unsigned int ncomp>
  class HostIrregData<BOUNDARY, data_t, ncomp>: public BaseHostIrregData<BOUNDARY, data_t, ncomp> 
  {
  public:

    typedef EBIndex<BOUNDARY> index_t;
    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    HostIrregData()
    {
    }

    HostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {
    }

    virtual bool copyContains(const index_t& a_ind, const Box & a_box) const 
    {

      bool retval = false;
      if(a_box.contains(a_ind.m_pt))
        {
          retval = this->contains(a_ind);
        }
      return retval;
    }

    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector<EBIndex<BOUNDARY> > >indicies(new vector< EBIndex<BOUNDARY> >() );
      *indicies = a_graph.getCutFaces(interBox);
      BaseHostIrregData<BOUNDARY, data_t, ncomp>::define(indicies);

    }
  };


  ///
  template <class data_t,  unsigned int ncomp>
  class HostIrregData<XFACE, data_t, ncomp>: public BaseHostIrregData<XFACE, data_t, ncomp> 
  {
  public:

    typedef EBIndex<XFACE> index_t;
    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    HostIrregData()
    {
    }

    virtual bool copyContains(const index_t& a_ind, const Box & a_box) const
    {

      bool retval = false;
      if(a_box.contains(a_ind.m_pt) || a_box.contains(a_ind.m_pt - Point::Basis(0)))
        {
          retval = this->contains(a_ind);
        }
      return retval;
    }

    HostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector<EBIndex<XFACE> > > indicies(new vector< EBIndex<XFACE> >() );
      *indicies = a_graph.getIrregXFaces(interBox);
      BaseHostIrregData<XFACE, data_t, ncomp>::define(indicies);
    }
  };


  ///
  template <class data_t,  unsigned int ncomp>
  class HostIrregData<YFACE, data_t, ncomp>: public BaseHostIrregData<YFACE, data_t, ncomp> 
  {
  public:

    typedef EBIndex<YFACE> index_t;
    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    HostIrregData()
    {
    }

    virtual bool copyContains(const index_t& a_ind, const Box & a_box) const
    {

      bool retval = false;
      if(a_box.contains(a_ind.m_pt) || a_box.contains(a_ind.m_pt - Point::Basis(1)))
        {
          retval = this->contains(a_ind);
        }
      return retval;
    }

    HostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector<EBIndex<YFACE> > > indicies(new vector< EBIndex<YFACE> >() );
      *indicies = a_graph.getIrregYFaces(interBox);
      BaseHostIrregData<YFACE, data_t, ncomp>::define(indicies);
    }
  };

  ///
  template <class data_t,  unsigned int ncomp>
  class HostIrregData<ZFACE, data_t, ncomp>: public BaseHostIrregData<ZFACE, data_t, ncomp> 
  {
  public:

    typedef EBIndex<ZFACE> index_t;
    ///
    HostIrregData( const vector<index_t> & a_indicies)
    {
      define(a_indicies);
    }

    HostIrregData()
    {
    }

    virtual bool copyContains(const index_t& a_ind, const Box & a_box) const
    {

      bool retval = false;
      if(a_box.contains(a_ind.m_pt) || a_box.contains(a_ind.m_pt - Point::Basis(2)))
        {
          retval = this->contains(a_ind);
        }
      return retval;
    }

    HostIrregData( const Box& a_box, bool a_stackAllocation=true)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
#if DIM==3
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<ZFACE> > > indicies(new vector< EBIndex<ZFACE> >() );
      *indicies = a_graph.getIrregZFaces(interBox);
      BaseHostIrregData<ZFACE, data_t, ncomp>::define(indicies);
#endif
    }
  };


}

#endif
