
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
  template <CENTERING srcCenter, CENTERING dstCenter>
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
    EBStencil(string                     a_stencilName,
              string                     a_domainBCName[2*DIM],
              string                     a_ebbcName,
              const DisjointBoxLayout  & a_grids
              const Point              & a_ghost);


    template<unsigned int Csrc, unsigned int Cdst>
    void apply(EBBoxData<dstCenter, double, 1> >       &   a_src,
               const EBBoxData<srcCenter, double, 1> > &   a_dst,
               bool            a_initToZero = false,
               const double    a_scale = 1) const;
               
              
  protected:
    
    ///these 
    Stencil<double> m_regStencilInterior;
    Stencil<double> m_regStencilBCS[2*DIM];
    
    AggStencil< EBData<srcCenter double, 1>, EBData<dstCenter, double, 1> > m_irregSten;


    private:
    EBStencil();
    EBStencil(     const EBStencil<srcCenter, dstCenter>& a_inp);
    void operator=(const EBStencil<srcCenter, dstCenter>& a_inp);


  };
}
#include "Proto_EBStencilImplem.H"
#endif
