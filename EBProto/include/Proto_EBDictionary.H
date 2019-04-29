
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
#include "EBProto.H"

using std::shared_ptr;
using std::vector;

namespace Proto
{
  ///
  /**
     makes geometric description from an implicit function.
  */
  template <int order, typename T, CENTERING srcCenter, CENTERING dstCenter>
  class EBDictionary
  {

  public:
    typedef EBStencil<T, srcCenter, dstCenter>  ebstencil_t;

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
    EBDictionary(const shared_ptr<GeometryService<order> >  & a_geoserv,
                 const DisjointBoxLayout                    & a_grids,
                 Point  a_srcGhost,  Point  a_dstGhost);


    inline void registerStencil(string a_stencilName,
                                string a_domainBCName[2*DIM],
                                string a_ebbcName);


    ///get the stencil for a particular box at a particular level of refinement
    inline shared_ptr< ebstencil_t >
    getEBStencil(const string & a_stencilName, 
                 const int    & a_boxid);

    ///
    inline int numLevels() const 
    {
      return m_registerStencils.size();
    }

  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
    void  close()
    {
      m_geoserv = shared_ptr<GeometryService<order> >();
    }
  protected:


    struct stencilKey
    {
      string m_ebbcName;
      string m_stencilName;
    };
 
    /// one per box in the grid
    vector< map< shared_ptr< ebstencil_t, stencilKey> >    m_registeredStencils;
    
    shared_ptr<GeometryService<order> > m_geoserv;
    Point                               m_srcGhost;
    Point                               m_dstGhost;
    DisjointBoxLayout                   m_grids;
    bool                                m_isClosed;
  private:
    EBDictionary();
    EBDictionary(  const EBDictionary<order>& a_inp);
    void operator=(const EBDictionary<order>& a_inp);


  };
}
#include "Proto_EBDictionaryImplem.H"
#endif
