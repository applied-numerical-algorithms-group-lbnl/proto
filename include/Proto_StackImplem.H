#include <cassert>
#include <cstdint>
#include <new>


void
StackAlloc::
StackAlloc()
{
#ifdef PROTO_CUDA
  cudaMalloc(&m_stack, STACKALLOC);
#else
  m_stack     = malloc(STACKALLOC);
#endif
  *m_stack = 0;
  m_top = m_stack;
}

void
StackAlloc::
~StackAlloc()
{
#ifdef PROTO_CUDA
  cudaFree(m_stack);
#else
  free(m_stack);
#endif
}

void* StackAlloc::alloc(size_t a_sz)
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

void StackAlloc::free(void* a_pt) // probably best if this seg-faults on nullptr, indicates a logic error
{
  //remember how the increment was stored in the extra bit at the end of the allocation?
  //this is only technically correct if objects are deallocated in the reverse order of allocation
  m_top  -= *m_top;

  //not sure if i want this check
//  PR_assert(m_top >= m_stack);
}

