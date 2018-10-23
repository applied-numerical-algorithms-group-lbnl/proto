#ifndef _STACK_ALLOC_H_
#define _STACK_ALLOC_H_

#include <cstddef>

///  class Stack implements a memory allocation as a static stack array.  
/** Stack is meant to be the default memory allocation for automatic variables that have runtime storage requirements, 
like FArrayBox's dataPtr.  It relies on the automatic variable semantics of LIFO.  This is literally a static storage array and
a stack-based LIFO control system.  When code is built with OpenMP, the static arrays are threadprivate.
*/

#define STACKALLOC 1073741824
class StackAlloc
{
    /// base class constructor
public:
  StackAlloc();

  virtual ~StackAlloc();

  /**
     Allocate a dynamic memory arena of size a_sz.
     A pointer to this memory should be returned.
  */
  static void* alloc(size_t a_sz);

  /**
     virtual function for deleting the arena pointed to by a_pt.
  */
  static void free(void* a_pt);


private:
  static size_t stack[STACKALLOC/sizeof(size_t)]; // 1GB static array
  static size_t* top;
#pragma omp threadprivate(stack, top)
};
#include "StackImplem.H"
#endif //end include guard
