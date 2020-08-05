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

#ifdef PR_MPI
#include "mpi.h"
#endif
// SPMD functions.
using namespace std;
namespace Proto
{
#ifdef PR_MPI
  template <typename T=void>
  struct Proto_MPI
  {
    static MPI_Comm comm;
  };  
#else
  // this can be changed for debugging parallel code in serial
  extern int num_procs;
#endif
  template <typename T>
  MPI_Comm Proto_MPI<T>::comm = MPI_COMM_WORLD;
#define PR_MAX_MPI_MESSAGE_SIZE 30*1024*1024  
  extern unsigned long long PR_MaxMPISendSize;
  extern unsigned long long PR_MaxMPIRecvSize;
#if 0
  inline int reportMPIStats()
  {
    pout()<<"Proto message size limit:"<< PR_MAX_MPI_MESSAGE_SIZE<<"\n"
          <<"Max send message size:"<<PR_MaxMPISendSize<<"\n"
          <<"Max recv message size:"<<PR_MaxMPIRecvSize<<std::endl;
  return 0;
  };
#endif
/// local process ID
/**
    Returns the ID of the locally running process in the range 0 <=
    procID() < numProc().  This has no relation to the operating system
    pid.  There is always a procID() == 0.  */
  inline int procID()
  {
    int retval;
    MPI_Comm_rank(Proto_MPI<void>::comm, &retval);
    return retval;
  };

  inline int PRprocID()
  {
    return procID();
  }
  
/// number of parallel processes
/**
   Returns the number of parallel processes running.
    Always returns at least 1. */
  inline unsigned int numProc()
  {
    static int ret = -1;
    if (ret == -1)
      {
        MPI_Comm_size(Proto_MPI<void>::comm, &ret);
      }
    return ret;
  };
///all ranks wait here to sync-up
/**
   All MPI ranks wait here to sync-up.  Calls MPI_Barrier(comm)  */
  inline void barrier()
  {
    MPI_Barrier(Proto_MPI<>::comm);
  };

}
#endif
