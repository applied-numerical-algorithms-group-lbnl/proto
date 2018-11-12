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
    
    Point dirlo = -Point::Basis(a_idir, a_nghost);
    Point dirhi =  Point::Basis(a_idir, a_nghost);

    Point unitdirlo = -Point::Basis(a_idir, 1);
    Point unitdirhi =  Point::Basis(a_idir, 1);
    Box dstbxlo = domain.edge(dirlo);
    Box dstbxhi = domain.edge(dirhi);

    //these are swapped because you have to move in different 
    //directions to get the periodic image
    Point shifthi = unitdirlo*domain.size(a_idir);
    Point shiftlo = unitdirhi*domain.size(a_idir);

    Box srcbxlo = dstbxlo.shift(shiftlo);
    Box srcbxhi = dstbxhi.shift(shifthi);

    a_phi.copy(a_phi, srcbxlo, 0, dstbxlo, 0, C);
    a_phi.copy(a_phi, srcbxhi, 0, dstbxhi, 0, C);
  }


} //end Proto namespace
#endif //end include guard
