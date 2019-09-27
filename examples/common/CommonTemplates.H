#ifndef _CommonTemplates_H_
#define _CommonTemplates_H_
#include "Proto.H"

namespace protocommon
{
  ///enforce periodic boundary conditions as if the box of the BoxData is the domain grown by nghost
  template <class T=double, unsigned int C>
  inline void
  enforceSGBoundaryConditions(BoxData<T, C>& a_phi,
                              int            a_nghost,
                              int            a_idir)
  {
    Box domain = a_phi.box().grow(-a_nghost);
    
    Point domlo = domain.low();
    Point domhi = domain.high();

    int domsize = domain.size(a_idir);
    {
      //high copies to low
      Box srcBoxStart = domain.adjCellHi(a_idir, a_nghost);
      Box srcBox = srcBoxStart.shift(-Point::Basis(a_idir, a_nghost));

      Point destShift = -Point::Basis(a_idir, domsize);
      a_phi.copyTo(a_phi, srcBox, destShift);
    }
    {
      //low copies to high
      Box srcBoxStart = domain.adjCellLo(a_idir, a_nghost);
      Box srcBox = srcBoxStart.shift(Point::Basis(a_idir, a_nghost));

      Point destShift = Point::Basis(a_idir, domsize);
      a_phi.copyTo(a_phi, srcBox, destShift);
    }
  }


} //end Proto namespace
#endif //end include guard
