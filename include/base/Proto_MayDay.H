#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _Proto_MAYDAY_H_
#define _Proto_MAYDAY_H_
#include "Proto_parstream.H"

#define CH_DEFAULT_ERROR_CODE 255
namespace Proto
{
/// Error-reporting Functions
/**
   This class consists entirely of static error-reporting functions that are
   used throughout the MayDay implementation together with a version string.
   In the best of all possible worlds we'd use a namespace, but we aren't
   there, yet :-) Note that all the output functions add a terminating
   exclamation mark, so there is no need to add any punction to sentences
   that get output explicitly.
*/
  template <typename T=void>
  class MayDay
  {
  public:
    
    static bool s_debugSpew;
    /// Print out message to cerr and exit with the specified exit code.
    [[noreturn]] inline static void Error(const char * const a_msg = m_nullString, int m_exitCode = CH_DEFAULT_ERROR_CODE);
    
    /// Print out message to cerr and exit via abort() (if serial) or MPI_Abort() (if parallel).
  /** The exit code depends on the compiler and MPI implementation.
   */
  inline static void Abort(const char * const a_msg = m_nullString);

  /// Print out message to cerr and continue
  inline static void Warning(const char * const a_msg = m_nullString);

  static const char * const version;

protected:
  /// A null string used as a default argument by the other member functions
  static const char * const m_nullString;
};
#include "implem/Proto_MayDayImplem.H"
} // end Proto namespace.
#endif
