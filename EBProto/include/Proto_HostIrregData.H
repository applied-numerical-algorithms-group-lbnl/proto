
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
    BaseHostIrregData( const Box& a_box)
    {

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
    inline std::size_t size() const {return m_indicies->size()*ncomp;};

    ///
    shared_ptr< vector<index_t> > getIndicies() const
    {
      return m_indicies;
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

    shared_ptr<vector<index_t> >  m_indicies;
    shared_ptr<data_t>            m_data;
    data_t*                       m_rawPtr;
  };


  ///
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class HostIrregData: public BaseHostIrregData<cent, data_t, ncomp> 
  {
  public:

    HostIrregData( const Box& a_box)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph);
  };



  ///
  template <class data_t,  unsigned int ncomp>
  class HostIrregData<CELL, data_t, ncomp>: public BaseHostIrregData<CELL, data_t, ncomp> 
  {
  public:

    HostIrregData( const Box& a_box)
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

    HostIrregData( const Box& a_box)
    {
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

    HostIrregData( const Box& a_box)
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

    HostIrregData( const Box& a_box)
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

    HostIrregData( const Box& a_box)
    {
    }
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<ZFACE> > > indicies(new vector< EBIndex<ZFACE> >() );
      *indicies = a_graph.getIrregYFaces(interBox);
      BaseHostIrregData<ZFACE, data_t, ncomp>::define(indicies);
    }
  };


}

#endif
