#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _IFSLICERIMPLEM_H_
#define _IFSLICERIMPLEM_H_

#include <iostream>

#include "NamespaceHeader.H"

// empty constructor
template <int dim> IFSlicer<dim>::IFSlicer()
{
  //set to undefined values
  m_IFSlicer     = NULL;
  m_fixedComp  = LARGEINTVAL;
  m_fixedValue = LARGEREALVAL;
}

//copy constructor
template <int dim >IFSlicer<dim>::IFSlicer(const IFSlicer<dim> & a_IFSlicer)
{
  setParams(a_IFSlicer.m_IFSlicer,
            a_IFSlicer.m_fixedComp,
            a_IFSlicer.m_fixedValue);
}

template <int dim> IFSlicer<dim>::IFSlicer(IFSlicer<dim+1> const * a_IFSlicer,
                                           const int             & a_fixedComp,
                                           const Real            & a_fixedValue)
{
  setParams(a_IFSlicer,
            a_fixedComp,
            a_fixedValue);
}

template <int dim> void IFSlicer<dim>::setParams(IFSlicer<dim+1> const * a_IFSlicer,
                                                 const int             & a_fixedComp,
                                                 const Real            & a_fixedValue)
{
  // copy the inputs
  m_IFSlicer   = a_IFSlicer;
  m_fixedComp  = a_fixedComp;
  m_fixedValue = a_fixedValue;
}

// Destructor
template <int dim> IFSlicer<dim>::~IFSlicer()
{
  //lower dimension doesn't own the objects referred to by any pointer.
}

template<int dim >Real IFSlicer<dim>::value(const IndexTM<int, dim> & a_partialDerivative,
                                            const IndexTM<Real,dim> & a_point) const
{
  Real retval = LARGEREALVAL;

  //use fixedComp and fixedValue to pull partialDerivative and point into one higher dimension
  IndexTM<int, dim + 1> partialDerivative;
  IndexTM<Real,dim + 1> point;
  for (int idir = 0; idir < dim +1; ++idir)
    {
      if (idir < m_fixedComp)
        {
          partialDerivative[idir] = a_partialDerivative[idir];
          point            [idir] = a_point            [idir];
        }
      else if (idir > m_fixedComp)
        {
          partialDerivative[idir] = a_partialDerivative[idir - 1];
          point            [idir] = a_point            [idir - 1];
        }
      else
        {
          partialDerivative[m_fixedComp] = 0;
          point            [m_fixedComp] = m_fixedValue;
        }
    }

  //evaluate at one higher dimension
  retval = m_IFSlicer->value(partialDerivative,point);

  return retval;
}

template <int dim> void IFSlicer<dim>::print(ostream& a_out) const
{
  a_out << "fixed component  = " << m_fixedComp  << "\n";
  a_out << "fixed value      = " << m_fixedValue << "\n";
  a_out << "At one dimension higher: " << dim + 1 << "\n";
  this->m_IFSlicer.print(a_out);
}

template <int dim> ostream& operator<<(ostream             & a_out,
                                       const IFSlicer<dim> & a_IFSlicer)
{
  a_IFSlicer.print(a_out);
  return a_out;
}
#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IFSlicer.H"

#include "NamespaceHeader.H"

// empty constructor (dim == GLOBALDIM)
IFSlicer<GLOBALDIM>::IFSlicer()
{
  m_implicitFunction = NULL;
}

//copy constructor
IFSlicer<GLOBALDIM>::IFSlicer(const IFSlicer<GLOBALDIM> & a_IFSlicer)
{
  m_implicitFunction = a_IFSlicer.m_implicitFunction->newImplicitFunction();
}

// constructor (dim == GLOBALDIM)
IFSlicer<GLOBALDIM>::IFSlicer(const BaseIF & a_implicitFunction)
{
  m_implicitFunction = a_implicitFunction.newImplicitFunction();
}

// Destructor (dim == GLOBALDIM)
IFSlicer<GLOBALDIM>::~IFSlicer()
{
  if (m_implicitFunction != NULL)
    {
      delete m_implicitFunction;
    }
}

Real IFSlicer<GLOBALDIM>::value(const IndexTM<int, GLOBALDIM> & a_partialDerivative,
                                const IndexTM<Real,GLOBALDIM> & a_point) const
{
  return m_implicitFunction->value(a_partialDerivative,a_point);
}

void IFSlicer<GLOBALDIM>::print(ostream& a_out) const
{
  MayDay::Abort("Not implemented");
  // m_implicitFunction->print(a_out);
}

#include "NamespaceFooter.H"

#include "NamespaceFooter.H"

#endif
