
#ifndef _Proto_EBStencil_H_
#define _Proto_EBStencil_H_


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


using std::shared_ptr;
using std::vector;

namespace Proto
{
  ///
  /**
  */
  template <class srcIndex, class dstIndex>
  class EBStencil
  {

  public:

    ///
    /**
       regular stencil gets applied in the interior of the domain.
       The aggstencil overwrites it near the EB.
       The domainBC stencils fill ghost cells outside the domain.  
       The vector length needs to be the same as the one-sided span of the stencil.
     */
    EBStencil(const Stencil<double>                             & a_regularStencil,
              const vector<Stencil<double> >                     a_domainBCStencils[2*DIM],
              const vector<shared_ptr<AggStencil<srcIndex, dstIndex> > >& a_irregStencil,
              const DisjointBoxLayout                           & a_grids,
              const Box                                         & a_regApplyBox);


    void apply
              
  protected:
    
    ///these 
    vector<shared_ptr<LevelData<VoluData> >   m_voludata;
    vector<shared_ptr<LevelData<FaceData> >   m_facedata[2*DIM]

  private:
    EBStencil();
    EBStencil(     const EBStencil<srcIndex, dstIndex>& a_inp);
    void operator=(const EBStencil<srcIndex, dstIndex>& a_inp);


  };
}
#include "Proto_GeometryServiceImplem.H"
#endif
