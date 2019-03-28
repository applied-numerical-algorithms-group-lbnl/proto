#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NOTATION_H_
#define _NOTATION_H_

#include <map>

#include "IndexTM.H"
#include "SPACE.H"

#include "NamespaceHeader.H"




// General constants
#define LARGEINTVAL -999999999
#define LARGEREALVAL 123456789.0
#define TOLERANCE            1.0e-10
#define MACHINEPRECISION     1.0e-14
#define REFINEMENTTHRESHOLD  1.0e-14

// These terms are possible values for cornerSigns
#define OUT -1
#define ON   0
#define IN   1

// These are for indexing into boundary maps
#define BDID_DIR   0
#define BDID_HILO  1

// Ordering for map
template<typename T> struct LexLT
{
  bool operator()(const T& a_i1,
                  const T& a_i2) const
  {
    return a_i1.lexLT(a_i2);
  }
};

enum EBorVol
{
  EBMoment  = 0,
  VolMoment = 1
};

// N-tuples of integers and reals
typedef IndexTM<int,2>            Iv2;
typedef IndexTM<int,GLOBALDIM>    IvgDim;
typedef IndexTM<int,GLOBALDIM-1 > IvgLess1;
typedef IndexTM<Real,GLOBALDIM>   RvgDim;

#include "NamespaceFooter.H"

#endif
