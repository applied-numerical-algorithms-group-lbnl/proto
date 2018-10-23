#include <cassert>
#include <cstdint>
#include <new>

size_t  StackAlloc::stack[STACKALLOC/sizeof(size_t)] ={0}; // 2GB static array, first entry set to 0
size_t* StackAlloc::top = StackAlloc::stack;


void* StackAlloc::alloc(size_t a_sz)
{
  void* rtn = (void*)(top+1);
  assert(a_sz>0);
  size_t increment = (2 + (a_sz-1)/sizeof(size_t))*sizeof(size_t); // safe round up divide and move ahead to next 
  top+= increment;
  if((uintptr_t)top > ((uintptr_t)stack+STACKALLOC)) throw std::bad_alloc();
  *top=increment;
  //std::cout<<"a_sz, increment, ptr: "<<a_sz<<" "<<increment<<" "<<(size_t)rtn<<"\n";
  return rtn;
}

void StackAlloc::free(void* a_pt) // probably best if this seg-faults on nullptr, indicates a logic error
{
  top -= *top;
}

