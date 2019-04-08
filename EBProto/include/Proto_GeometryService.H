
#ifndef _Proto_WRAPPEDGSHOP_H_
#define _Proto_WRAPPEDGSHOP_H_


#include "Proto_RealVect.H"
#include "Proto_Box.H"
#include "Proto_Point.H"
#include "Proto_LSquares.H"
#include "Proto_BaseIF.H"
#include "Proto_IrregNode.H"
#include "Proto_CutCellMoments.H"
#include "Proto_HostBoxData.H"
#include "Proto_PointSet.H"
using std::shared_ptr;
using std::vector;

namespace Proto
{

  ///
  /**
     makes geometric description from an implicit function.
  */
  template <int order>
  class GeometryService
  {

  public:
    typedef IrregNode<order> inode;
    ///
    /**
       This class will refine a cell if it is between min and max refinements.
       If it violates bounds, that triggers refinement until max refinement.
    */
    GeometryService(const shared_ptr<BaseIF>  & a_baseIF,
                    const RealVect            & a_origin,
                    const double              & a_dx,
                    const Box                 & a_domain);


    ///
    /**
       Define the internals of the input ebisRegion.
    */
    inline void fillGraph(HostBoxData<int>          &   a_regIrregCovered,
                          vector<inode>             &   a_nodes,
                          const Box                 &   a_validRegion,
                          const Box                 &   a_ghostRegion) const;


  protected:

    ///
    ~GeometryService()
    {
    }

    
    bool onBoxBoundary(const Point        & a_iv, 
                       const Box            & a_box,
                       const int            & a_dir,
                       const Side::LoHiSide & a_sd) const;



    //converts a RealVect in physical coordinates to a RealVect in coordinates relative to a cell center
    inline RealVect convert2RelativeCoord(const RealVect& a_rVect)const;

    //converts a IndexTM<double,DIM> in physical coordinates to a
    //RealVect in coordinates relative to a cell center
    inline RealVect convert2RelativeCoord(const IndexTM<double,DIM>& a_rVect)const;

    inline void
    fillNewNode(inode                    &     a_node,
                const PointSet           &     a_ivsIrreg,
                const Point              &     a_iv) const;


  private:
    inline void
    fixRegularCellsNextToCovered(vector<inode>       & a_nodes, 
                                 HostBoxData<int>    & a_regIrregCovered,
                                 const Box           & a_validRegion,
                                 const Point         & a_iv) const;
    inline void
    getFullNodeWithCoveredFace(inode                   & a_newNode, 
                               const HostBoxData<int>  & a_regIrregCovered,
                               const Point             & a_iv);

  
  
    double                m_threshold;

    //origin
    RealVect            m_origin;

    int                 m_order;
    int                 m_degreeP;
    double              m_dx;

    Box m_domain;

    shared_ptr<BaseIF>       m_baseIF;
  

    inline void fillArc(vector<int>                          a_arc[DIM],
                        CutCellMoments<DIM>       &     a_cutCellMoments,
                        const int                      &     a_hilo,
                        const PointSet                 &     a_ivsIrreg,
                        const Point                  &     a_curriv) const;


    //stuff disallowed for all the usual reasons.
  private:
    GeometryService();

    GeometryService(const GeometryService& a_workshopin);

    void operator=(const GeometryService& a_workshopin);


  };
}
#include "Proto_GeometryServiceImplem.H"
#endif
