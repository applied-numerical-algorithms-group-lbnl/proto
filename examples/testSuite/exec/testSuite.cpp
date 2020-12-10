#include <iostream>
#include <test_forall.cu>
#include <test_reduction.cu>
#include <test_boxdata_operators.cu>
//#include <test_stack.cu>
#ifdef PROTO_CUDA
#include <test_fusion_bc.cu>
#endif


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
#ifdef PROTO_CUDA
  cudaSetDevice(4);
#endif

#ifdef PROTO_CUDA
  do_test("test_fusion_bc",   run_test_fusion_bc); 
#endif
  do_test("test_forall",      run_test_forall); 
  do_test("test_forall_i",    run_test_forall_p); 
  do_test("test_forall_p",    run_test_forall_i); 
//  do_test("test_stack_using", run_test_stack_using); 
//  do_test("test_stack_free",  run_test_stack_free); 
//  do_test("test_stack_empty",  run_test_stack_empty); 
//  do_test("test_stack_reset",  run_test_stack_empty); 
  do_test("test_reduction_min_linear_init_1",  test_reduction_min_linear_init_1); 
  do_test("test_reduction_min_linear_init_minus_2",  test_reduction_min_linear_init_minus_2); 
  do_test("test_reduction_max_linear_init_1",  test_reduction_max_linear_init_1); 
  do_test("test_reduction_max_linear_init_minus_2",  test_reduction_max_linear_init_minus_2); 
  do_test("test_reduction_abs_linear_init_1",  test_reduction_abs_linear_init_1); 
  do_test("test_reduction_abs_linear_init_minus_2",  test_reduction_abs_linear_init_minus_2); 
  do_test("test_reduction_reset_min", test_reduction_reset_min); 
  do_test("test_reduction_reset_max", test_reduction_reset_max); 
  do_test("test_reduction_reset_abs", test_reduction_reset_abs); 
  do_test("test_boxdata_operators_set_value_zero",test_boxdata_operators_set_value_zero);
  do_test("test_boxdata_operators_set_value_two",test_boxdata_operators_set_value_two);
  do_test("test_boxdata_operators_set_value_minus_three",test_boxdata_operators_set_value_minus_three);
  do_test("test_boxdata_operators_add_one_and_two",test_boxdata_operators_add_one_and_two);
  do_test("test_boxdata_operators_subtract_one_and_two",test_boxdata_operators_subtract_one_and_two);
  do_test("test_boxdata_operators_multiply_one_and_two",test_boxdata_operators_multiply_one_and_two);
  do_test("test_boxdata_operators_divide_one_and_two",test_boxdata_operators_divide_one_and_two);
  do_test("test_boxdata_operators_add_neg",test_boxdata_operators_add_neg);
  do_test("test_boxdata_operators_subtract_neg",test_boxdata_operators_subtract_neg);
  do_test("test_boxdata_operators_multiply_neg",test_boxdata_operators_multiply_neg);
  do_test("test_boxdata_operators_divide_neg",test_boxdata_operators_divide_neg);
  do_test("test_boxdata_operators_copy_full",test_boxdata_operators_copy_full);
  do_test("test_boxdata_operators_copy_to_smaller_box_data",test_boxdata_operators_copy_to_smaller_box_data);
  do_test("test_boxdata_operators_copy_to_smaller_box",test_boxdata_operators_copy_to_smaller_box);
  return 0;  
}
