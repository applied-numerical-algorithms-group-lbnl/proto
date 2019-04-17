
#ifndef _Proto_EBDictionary_H_
#define _Proto_EBDictionary_H_


#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "Proto_IrregNode.H"
#include "Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "Proto_HostIrregData.H"
#include "Proto_PointSet.H"
#include "Proto_LevelData.H"
#include "Proto_EBGraph.H"
#include "Proto_EBStencil.H"

using std::shared_ptr;
using std::vector;

namespace Proto
{
  ///
  /**
     makes geometric description from an implicit function.
  */
  template <int order>
  class EBDictionary
  {

  public:

    typedef IndexedMoments<DIM  , order> IndMomDIM;
    typedef IndexedMoments<DIM-1, order> IndMomSDMinOne;
    typedef MomentIterator<DIM  , order> MomItDIM;
    typedef MomentIterator<DIM-1, order> MomItSDMinOne;
    typedef IrregNode<order> inode;
    typedef IrregData<IndMomDIM     , VolIndex , 1> VoluData;
    typedef IrregData<IndMomSDMinOne, FaceIndex, 1> FaceData;

    ///
    /**
       a_ghost is the number of ghost cells your data will hold.
       It must be at least as big as the stencil span in each direction..
    */
    EBDictionary(const GeometryService<order> & a_geoserv,
                 const RealVect               & a_origin,
                 const double                 & a_dx,
                 const DisjointBoxLayout      & a_finestGrids,
                 const Point                  & a_ghost,
                 int  a_maxCoarsen = 0)


    inline void registerStencil(string a_stencilName);

    ///after this is done, you cannot register any more stencils.  This will delete the moment information.
    inline void close();


    inline shared_ptr<EBStencil> getEBStencil(const string & a_stencilName, 
                                              const Box    & a_domain);

    int numLevels() const 
    {
      return m_registerStencils.size();
    }
  protected:

    ///these vectors are of the length of number of coarsenings
    vector<shared_ptr<LevelData<VoluData> > >       m_voludata;

    ///
    vector<shared_ptr<LevelData<FaceData> > >       m_facedata[2*DIM];

    ///
    vector<shared_ptr<LevelData<EBGraph> >  m_graphs;

    ///
    vector< map<shared_ptr<EBStencil>, string> >   m_registeredStencils;

  private:
    EBDictionary();
    EBDictionary(  const EBDictionary<order>& a_inp);
    void operator=(const EBDictionary<order>& a_inp);


  };
}
#include "Proto_EBDictionaryImplem.H"
#endif
