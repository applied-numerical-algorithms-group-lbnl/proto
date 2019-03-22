#ifndef _proto_BASEIF_H_
#define _proto_BASEIF_H_

#include "Proto_RealVect.H"
#include "Proto_IndexTM.H"
#include "Proto_Point.H"
#include "Proto_IndexedMoments.H"


class Proto
{
  ///
  /**
     This is the base class for an implicit function specification of geometry.
     All that is needed is a constructor/destructor, a method to give the value
     of the function at any point in space (in 2D or 3D), and a factory method
     (these will probably all be very similar).
  */
  class BaseIF
  {

  public:

    /// Default constructor
    BaseIF()
    {
    }

    /// Default destructor
    virtual ~BaseIF()
    {
    }

    ///
    /**
       Return the value of the function at a_point.  When delineating a domain,
       the level set value=0 represents the boundary and value<0 is inside the
       fluid.
    */
    virtual Real value(const RealVect& a_point) const = 0;

    ///return the partial derivative at the point
    virtual Real derivative(const  Point  & a_deriv,
                            const RealVect& a_point) const = 0;

    virtual Real value(const IndexTM<int,DIM> & a_partialDerivative,
                       const IndexTM<Real,DIM>& a_point) const
    {
      Real retval= 0;
      if (a_partialDerivative == IndexTM<int,DIM>::Zero)
      {
        retval = value(a_point);
      }
      else
      {
        retval = derivative(a_partialDerivative, a_point);
      }
      return retval;
    }

    
    ///
    /**
       Return a newly allocated derived class.  The responsibility
       for deleting the memory is left to the calling function.
    */
    virtual BaseIF* newImplicitFunction() const = 0;

  };

}
#endif
