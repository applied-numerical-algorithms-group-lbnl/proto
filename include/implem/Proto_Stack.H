#ifndef _STACK_ALLOC_H_
#define _STACK_ALLOC_H_

#include <cstddef>
#include <stack>
#include <iostream>
#include <memory>
using std::shared_ptr;
///  class Stack implements a memory allocation as chunks of a big array.
/** Stack is meant to be the default memory allocation for automatic variables that have runtime storage requirements, 
    like FArrayBox's dataPtr.  It relies on the automatic variable semantics of LIFO.  This is literally a static storage array and
    a stack-based LIFO control system.  When code is built with OpenMP, the static arrays are threadprivate.
*/

#define STACKALLOC 1073741824
//#define STACKALLOC 2147483648

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
    m_top = m_stack;
    m_maxusage= 0;
  }


  inline ~StackAlloc()
  {
#ifdef PROTO_CUDA
    cudaFree(m_stack);
#else
    free(m_stack);
#endif
//    std::cout << "max stack memory usage at destruction = " << m_maxusage << std::endl;
  }

  /**
     Allocate a dynamic memory arena of size a_sz.
     A pointer to this memory should be returned.
  */
  inline void* alloc(size_t a_sz)
  {
    void* rtn = (void*)(m_top);

    PR_assert(a_sz>0);
    //compute the safe increment into the stack and reset the top value
    //size_t increment = (2 + (a_sz-1)/sizeof(size_t))*sizeof(size_t); // safe round up divide and move ahead to next
     size_t increment = a_sz;
//    std::cout << "alloc increment = " << increment << std::endl;
    m_top += increment;
    m_increments.push(increment);

    size_t curusage = (m_top - m_stack);
    if(curusage > m_maxusage)
    {
      m_maxusage = curusage;
    }
    
    //make certain we have not gone off the end of the world
    //if((uintptr_t)m_top > ((uintptr_t)m_stack+STACKALLOC)) throw std::bad_alloc();
    if(curusage > STACKALLOC) throw std::bad_alloc();

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
    if(!m_increments.empty())
    {
//    size_t  decrement = m_increments.top();
//    std::cout << "free decrement = " << decrement << std::endl;
      m_top  -= m_increments.top();
      m_increments.pop();
    }
  }

  inline void reset()
  {
    while (!m_increments.empty())
    {
      m_increments.pop();
    }
    m_top = m_stack;
  }
private:
  size_t*            m_stack;
  size_t*            m_top;
  size_t             m_maxusage;
  std::stack<size_t> m_increments;

//#pragma omp threadprivate(m_stack, m_top, m_maxusage, m_increments)

};

#endif //end include guard
