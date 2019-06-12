#include <cassert>
#include <cstdint>
#include <new>


void* StackAlloc::alloc(size_t a_sz)
{
  shared_ptr<stackData> sdat = getStackData();
  static size_t* stack = sdat->stack;
  static size_t* top   = sdat->top;

  
  void* rtn = (void*)(top+1);

  PR_assert(a_sz>0);
  size_t increment = (2 + (a_sz-1)/sizeof(size_t))*sizeof(size_t); // safe round up divide and move ahead to next 
  top+= increment;
  if((uintptr_t)top > ((uintptr_t)stack+STACKALLOC)) throw std::bad_alloc();
  *top=increment;
  //std::cout<<"a_sz, increment, ptr: "<<a_sz<<" "<<increment<<" "<<(size_t)rtn<<"\n";

  sdat->top = top;
  sdat->sizes.push(a_sz);
  return rtn;
}

void StackAlloc::free(void* a_pt) // probably best if this seg-faults on nullptr, indicates a logic error
{
  
  shared_ptr<stackData> sdat = getStackData();
  if(!sdat->sizes.empty())
  {
    size_t latest = sdat->sizes.top();
    sdat->top -= latest;
    sdat->sizes.pop();
  }
}

