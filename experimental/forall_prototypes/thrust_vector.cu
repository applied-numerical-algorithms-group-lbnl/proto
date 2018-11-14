#include <iostream>
#include <cstdio>
#include <vector>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>


int main ( int argc, char** argv)
{
  int nsize = 16;
  std::vector<int> h_vec_std(nsize, 7);
  
  thrust::host_vector<int> h_vec_thr1(nsize);
  thrust::copy(thrust::host, h_vec_std.begin(), h_vec_std.end(), h_vec_thr1.begin());

  thrust::device_vector<int> d_vec(h_vec_thr1);
  thrust::host_vector<int> h_vec_thr = d_vec;

/**/

  for(int ivec = 0; ivec < nsize; ivec++)
  {
    std::cout << h_vec_thr[ivec] << " ";
  }
  std::cout << std::endl;
/**/
  return 0;
}

