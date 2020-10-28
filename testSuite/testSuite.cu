#include <iostream>
#include <test_fusion_bc.cu>
#include <test_forall.cu>
#include <test_stack.cu>


template<typename Func>
void do_test(std::string a_str, Func &fun)
{
  std::cout << "Do " << a_str << std::endl;
  bool b = fun();
  if(b) std::cout << "-> passed " << a_str << std::endl;
  else std::cout << "-> failed " << a_str << std::endl;
}

int main()
{
  cudaSetDevice(1);
  do_test("test_fusion_bc",   run_test_fusion_bc); 
  do_test("test_forall",      run_test_forall); 
  do_test("test_forall_i",    run_test_forall_p); 
  do_test("test_forall_p",    run_test_forall_i); 
  do_test("test_stack_using", run_test_stack_using); 
  do_test("test_stack_free",  run_test_stack_free); 
  return 0;  
}
