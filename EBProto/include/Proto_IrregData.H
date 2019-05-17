
#ifndef _Proto_IrregDATA_H_
#define _Proto_IrregDATA_H_

#include <cmath>
#include <cstdlib>
#include <memory>
#include "Proto_Box.H"
#include "EBProto.H"

using std::vector;
using std::shared_ptr;
namespace Proto
{


///
/**
*/
  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class BaseIrregData
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    BaseIrregData()
    {
    }

    ///
    virtual ~BaseIrregData()
    {
    }

    ///
    BaseIrregData( const shared_ptr< vector<index_t> >& a_indicies)
    {
      define(a_indicies);
    }

    ///needed for leveldata
    BaseIrregData(const Box& a_box)
    {  }

    ///
    inline void define( const shared_ptr< vector< index_t > >& a_indicies);

    ///
    inline std::size_t size() const {return m_indicies->size()*ncomp;};

    ///
    inline bool defined() const {return bool(m_data);};
   

    ///
    inline void setVal(const data_t& a_val);
  

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
      PROTO_ASSERT(found, "vof not found");
      ivec += a_ivar*m_indicies->size();
      return ivec;
    }



  protected:

    shared_ptr<vector<  index_t > >        m_indicies;
    shared_ptr<data_t>         m_data;
    data_t*                    m_rawPtr;

  };


  template <CENTERING cent, class data_t,  unsigned int ncomp>
  class IrregData: public BaseIrregData<cent, data_t, ncomp>
  {
  public:
    typedef EBIndex<cent> index_t;
    ///
    IrregData()
    {
    }

    virtual ~IrregData()
    {
    }

    ///
    IrregData( const shared_ptr< vector<index_t> >& a_indicies)
    {
      define(a_indicies);
    }


    inline void define(const Box& a_box,  const EBGraph& a_graph);
  };

  ///template specialization for cell based data
  template < class data_t,  unsigned int ncomp>
  class IrregData <CELL, data_t, ncomp>: public BaseIrregData<CELL, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box,  const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();
      shared_ptr<vector<EBIndex<CELL> > > indicies(new vector< EBIndex<CELL> >() );
      *indicies = a_graph.getIrregVoFs(interBox);

      BaseIrregData<CELL, data_t, ncomp>::define(indicies);
    }
  };


  ///template specialization for boundary data
  template < class data_t,  unsigned int ncomp>
  class IrregData <BOUNDARY, data_t, ncomp>: public BaseIrregData<BOUNDARY, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<BOUNDARY> > > indicies(new vector< EBIndex<BOUNDARY> >() );
      *indicies = a_graph.getCutFaces(interBox);
      BaseIrregData<BOUNDARY, data_t, ncomp>::define(indicies);

    }
  };

  ///template specialization for xface data
  template < class data_t,  unsigned int ncomp>
  class IrregData <XFACE, data_t, ncomp>: public BaseIrregData<XFACE, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<XFACE> > > indicies(new vector< EBIndex<XFACE> >() );
      *indicies = a_graph.getIrregXFaces(interBox);
      BaseIrregData<XFACE, data_t, ncomp>::define(indicies);

    }

  };

  ///template specialization for yface data
  template < class data_t,  unsigned int ncomp>
  class IrregData <YFACE, data_t, ncomp>: public BaseIrregData<YFACE, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<YFACE> > > indicies(new vector< EBIndex<YFACE> >() );
      *indicies = a_graph.getIrregYFaces(interBox);
      BaseIrregData<YFACE, data_t, ncomp>::define(indicies);

    }

  };

#if DIM==3
  ///template specialization for zface data
  template < class data_t,  unsigned int ncomp>
  class IrregData <ZFACE, data_t, ncomp>: public BaseIrregData<ZFACE, data_t, ncomp>
  {
  public:
    inline void define(const Box& a_box, const EBGraph& a_graph)
    {
      Box interBox = a_box & a_graph.getDomain();

      shared_ptr<vector< EBIndex<ZFACE> > > indicies(new vector< EBIndex<ZFACE> >() );
      *indicies = a_graph.getIrregZFaces(interBox);
      BaseIrregData<ZACE, data_t, ncomp>::define(indicies);

    }
  };
#endif
}


#include "Proto_IrregDataImplem.H"

#endif
