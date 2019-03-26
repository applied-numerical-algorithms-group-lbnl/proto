
#ifndef _Proto_WRAPPEDGSHOP_H_
#define _Proto_WRAPPEDGSHOP_H_


#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "Proto_IrregNode.H"
#include "Proto_CutCellMoments.H"
using std::shared_ptr;
using std::vector;
namespace Proto
{

///
/**
   makes geometric description from an implicit function.
*/
  class WrappedGShop
  {

  public:
    ///
    /**
       This class will refine a cell if it is between min and max refinements.
       If it violates bounds, that triggers refinement until max refinement.
    */
    WrappedGShop(const shared_ptr<BaseIF>  & a_baseIF,
                 const RealVect            & a_origin,
                 const double              & a_dx,
                 const Box                 & a_domain);


    ///
    /**
       Define the internals of the input ebisRegion.
    */
    inline void fillGraph(BaseFab<int>      &   a_regIrregCovered,
                          Vector<IrregNode> &   a_nodes,
                          const Box         &   a_validRegion,
                          const Box         &   a_ghostRegion,
                          const Box         &   a_domain,
                          const RealVect    &   a_origin,
                          const double      &   a_dx) const



  protected:

    ///
    ~WrappedGShop()
    {
    }


    bool onBoxBoundary(const Point        & a_iv, 
                       const Box            & a_box,
                       const int            & a_dir,
                       const Side::LoHiSide & a_sd) const;

    /**
     */
    inline void
    computeVoFInternals(IrregNode                &     a_node,
                        const IntVectSet         &     a_ivsIrreg,
                        const Box      &     a_domain,
                        const RealVect           &     a_origin,
                        const double               &     a_dx,
                        const IntVect            &     a_iv) const;



    //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
    inline RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

    //converts a IndexTM<double,SpaceDim> in physical coordinates to a
    //RealVect in coordinates relative to a cell center
    inline RealVect convert2RelativeCoord(const IndexTM<double,SpaceDim>& a_rVect)const;

    inline void
    fillNewNode(IrregNode                &     a_node,
                const IntVectSet         &     a_ivsIrreg,
                const Box      &     a_domain,
                const RealVect           &     a_origin,
                const double               &     a_dx,
                const IntVect            &     a_iv) const;


  private:
    inline void
    fixRegularCellsNextToCovered(Vector<IrregNode>   & a_nodes, 
                                 BaseFab<int>        & a_regIrregCovered,
                                 const Box           & a_validRegion,
                                 const Box           & a_domain,
                                 const IntVect       & a_iv,
                                 const double          & a_dx) const;
    inline void
    getFullNodeWithCoveredFace(IrregNode            & a_newNode, 
                               const BaseFab<int>   & a_regIrregCovered,
                               const IntVect        & a_iv,
                               const double           & a_dx,
                               const Box            & a_domain) const;

  
  
    double                m_threshold;

    //origin
    RealVect            m_origin;

    int                 m_order;
    int                 m_degreeP;


    Box m_domain;

    shared_ptr<BaseIF>       m_baseIF;
  

    inline void fillArc(Vector<int>                          a_arc[SpaceDim],
                        CutCellMoments<SpaceDim>       &     a_cutCellMoments,
                        const int                      &     a_hilo,
                        const IntVectSet               &     a_ivsIrreg,
                        const IntVect                  &     a_curriv) const;


    //stuff disallowed for all the usual reasons.
  private:
    WrappedGShop();

    WrappedGShop(const WrappedGShop& a_workshopin);

    void operator=(const WrappedGShop& a_workshopin);


  };
#include "Proto_WrappedGShopImplem.H"
}
#endif
