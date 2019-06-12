#ifndef _STACK_ALLOC_H_
#define _STACK_ALLOC_H_

#include <cstddef>
#include <stack>
using std::shared_ptr;
///  class Stack implements a memory allocation as chunks of a big array.
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
  void* alloc(size_t a_sz);

  /**
     virtual function for deleting the arena pointed to by a_pt.
     rewinds top by size of last allocaiton.
  */
  void free(void* a_pt);

  void reset()
  {
    m_top = m_stack;
  }
private:
  size_t*  m_stack;
  size_t*  m_top;
  static shared_ptr<stackData> getStackData()
  {
    static shared_ptr<stackData> stackDat;
    static bool alloced = false;
    if(!alloced)
    {
      stackDat = shared_ptr<stackData>(new stackData());
      stackDat->top = 0;
#ifdef PROTO_CUDA
      cudaMalloc(&stackDat->stack, STACKALLOC);
#else
      stack = malloc(STACKALLOC);
#endif
      alloced = true;
    }
    return stackDat;
  }
private:


};
#include "Proto_StackImplem.H"
#endif //end include guard
