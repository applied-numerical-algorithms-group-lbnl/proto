#include <iostream>
#include <test_fusion_bc.cu>
#include <test_forall.cu>
#include <test_stack.cu>
#include <test_reduction.cu>


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
  do_test("test_stack_empty",  run_test_stack_empty); 
  do_test("test_stack_reset",  run_test_stack_empty); 
  do_test("test_reduction_min_linear_init_1",  test_reduction_min_linear_init_1); 
  do_test("test_reduction_min_linear_init_minus_2",  test_reduction_min_linear_init_minus_2); 
  do_test("test_reduction_max_linear_init_1",  test_reduction_max_linear_init_1); 
  do_test("test_reduction_max_linear_init_minus_2",  test_reduction_max_linear_init_minus_2); 
  return 0;  
}
