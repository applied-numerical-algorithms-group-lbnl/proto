#ifndef _STACK_ALLOC_H_
#define _STACK_ALLOC_H_

#include <cstddef>
#include <stack>
#include <memory>
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
  inline StackAlloc()
  {
#ifdef PROTO_CUDA
    cudaMalloc(&m_stack, STACKALLOC);
#else
    m_stack     = (size_t*)(malloc(STACKALLOC));
#endif
    *m_stack = 0;
    m_top = m_stack;
  }


  inline ~StackAlloc()
  {
    size_t maxusage = m_top - m_stack;
#ifdef PROTO_CUDA
    cudaFree(m_stack);
#else
    free(m_stack);
#endif
  }

  /**
     Allocate a dynamic memory arena of size a_sz.
     A pointer to this memory should be returned.
  */
  inline void* alloc(size_t a_sz)
  {
    //give the current starting address +1 as the return value (because m_top is reserved for the increment size)
    void* rtn = (void*)(m_top+1);

    PR_assert(a_sz>0);
    //compute the safe increment into the stack and reset the top value
    size_t increment = (2 + (a_sz-1)/sizeof(size_t))*sizeof(size_t); // safe round up divide and move ahead to next 
    m_top += increment;

    //make certain we have not gone off the end of the world
    if((uintptr_t)m_top > ((uintptr_t)m_stack+STACKALLOC)) throw std::bad_alloc();

    //set the final bit of memory to the size of the increment
    *m_top=increment;

    return rtn;
  }

  /**
     virtual function for deleting the arena pointed to by a_pt.
     rewinds top by size of last allocaiton.
  */
  inline void free(void* a_pt)
  {
    //remember how the increment was stored in the extra bit at the end of the allocation?
    //this is only technically correct if objects are deallocated in the reverse order of allocation
    m_top  -= *m_top;

    //not sure if i want this check
//  PR_assert(m_top >= m_stack);
  }

  inline void reset()
  {
    m_top = m_stack;
  }
private:
  size_t*  m_stack;
  size_t*  m_top;



};

#endif //end include guard
