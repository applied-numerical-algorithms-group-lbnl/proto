#ifndef _IFSLICER_H_
#define _IFSLICER_H_

#include <map>
using std::vector;

#include "Proto_IndexTM.H"
#include "Proto_BaseIF.H"

using std::endl;

namespace Proto
{

  template <int dim> class IFSlicer
  {
  public:

    /// null constructor
    /**
       null constructor
    */
    IFSlicer();

    /// copy constructor
    /**
       copy constructor
    */
    IFSlicer(const IFSlicer<dim> & a_IFSlicer);

    /// makes a reduced IFSlicer
    /**
       makes a reduced IFSlicer
    */
    IFSlicer(IFSlicer<dim+1> const * a_slicer,
             const int             & a_fixedComp,
             const Real            & a_fixedValue);

    void setParams(IFSlicer<dim+1> const * a_slicer,
                   const int             & a_fixedComp,
                   const Real            & a_fixedValue);

    /// destructor
    /**
       destructor
    */
    virtual ~IFSlicer();

    /// Return the partial derivative evaluated at a_point.
    /**
       Represent the partial derivative as a list of integers.
       Return the partial derivative evaluated at a_point.
    */
    virtual Real value(const IndexTM<int, dim> & a_partialDerivative,
                       const IndexTM<Real,dim> & a_point)const;

  protected:
    // The higher dimensional IFSlicer being further sliced
    IFSlicer<dim + 1> const * m_IFSlicer;

    // The component direction being held fixed and its value
    int  m_fixedComp;
    Real m_fixedValue;

  private:
    //
    // These functions are made private to prevent use of the default
    // functions provided by the C++ compiler.
    //
    IFSlicer<dim>& operator=(const IFSlicer<dim>&);

  };

  template<> class IFSlicer<DIM>
  {
  public:
    /// null constructor
    /**
       null constructor
    */
    IFSlicer();

    /// copy constructor
    /**
       copy constructor
    */
    IFSlicer(const IFSlicer<DIM> & a_IFSlicer);

    /// if dim == DIM, construct from BaseIF
    /**
       if dim == DIM, construct from BaseIF
    */
    IFSlicer(const BaseIF & a_implicitFunction);

    /// destructor
    /**
       destructor
    */
    virtual ~IFSlicer();

    /// Return the partial derivative evaluated at a_point.
    /**
       Represent the partial derivative as a list of integers.
       Return the partial derivative evaluated at a_point.
    */
    virtual Real value(const IndexTM<int, DIM> & a_partialDerivative,
                       const IndexTM<Real,DIM> & a_point) const;

  protected:
    // The implicit function being sliced.  At DIM this is just a
    // container class to start the slicing process.
    BaseIF* m_implicitFunction;

  private:
    //
    // These functions are made private to prevent use of the default
    // functions provided by the C++ compiler.
    //
    IFSlicer<DIM>& operator=(const IFSlicer<DIM>&);
  };

}

#include "Proto_IFSlicerImplem.H"

#endif
