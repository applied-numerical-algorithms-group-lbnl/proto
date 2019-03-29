
#ifndef _Proto_IRREGNODE_H_
#define _Proto_IRREGNODE_H_


#include "Proto_Point.H"
#include "Proto_RealVect.H"
#include "Proto_IndexedMoments.H"

using std::vector;
namespace Proto
{
///
/**
   Node for construction of geometric information.
*/
  template <int order>
  class IrregNode
  {
  public:
    typedef IndexedMoments<SpaceDim  , order> IndMomSpaceDim;
    typedef IndexedMoments<SpaceDim-1, order> IndMomSDMinOne;

    ///
    IrregNode()
    {;    }

    ///
    ~IrregNode()
    {;    }

    ///return the index into the arc vectors
    static inline int index(int a_idir, Side::LoHiSide a_side) ;

    ///
    Point     m_cell;


    ///
    /**
       Indicies into a_nodes to show connectivity.
       If the arc is to  an irregular
       cell, the index is the unique index of the vof
       in the cell.
       For arcs to regular cells, the arc = -2
       If the arc is to the domain boundary, arc  = -1.
    */
    vector<int>      m_arc[2*DIM];


    ///for covered cells
    inline void setMomentsToZero();

    //for regular cells
    inline void setMomentsToRegular(const double& a_dx);

    ///data for irregular nodes

    /// the moment at  the VoF associated with the  monomial with the input exponents
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_VoF(x^p y^q z^r dV)
    **/
    IndMomDIM m_volumeMoments;


    /// the moment at  the irregular face associated with the  monomial with the input exponents 
    /**
       Given VoF variables x, y, z,   p = mono(0), q = mono(1), r = mono(2),
       returns integral_over_irregular_area((x^p y^q z^r) dA)
    **/
    IndMomDIM m_EBMoments;
  
    IndMomDIM m_normalPartialDeriv[DIM];
    ///face centered moments
    /// the moment at  the face associated with the  monomial with the input exponents 
    /**
       Given face variables x, y,    p = mono(0), q = mono(1)
       returns integral_over_face_area((x^p y^q ) dA)
    **/
    IndMomSDMinOne m_faceMoments[2*DIM];


  };

}


#endif
