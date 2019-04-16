
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
       This class will refine a cell if it is between min and max refinements.
       If it violates bounds, that triggers refinement until max refinement.
    */
    EBDictionary(const GeometryService<order> & a_baseIF,
                 const RealVect               & a_origin,
                 const double                 & a_dx,
                 const Box                    & a_domain,
                 const int                    & a_maxGridSize,
                 bool a_saveMoments);


    inline void registerStencil(string a_stencilName);

    inline shared_ptr<EBStencil> getEBStencil(const string & a_stencilName, 
                                              const Box    & a_domain
  protected:
    
    ///these 
    vector<shared_ptr<LevelData<VoluData> >   m_voludata;
    vector<shared_ptr<LevelData<FaceData> >   m_facedata[2*DIM]

  private:
    EBDictionary();
    EBDictionary(  const EBDictionary<order>& a_inp);
    void operator=(const EBDictionary<order>& a_inp);


  };
}
#include "Proto_GeometryServiceImplem.H"
#endif
