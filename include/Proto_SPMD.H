#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef _PROTO_SPMD_H_
#define _PROTO_SPMD_H_
#include <array>
#include "Proto.H"

// dummy SPMD functions for debugging.
using namespace std;
namespace Proto
{
  int nProcs()
  {return 1;};
  int myProc()
  {return 0;};
  int procID()
  {
    return 0;
  };
}
#endif
