

#include <stdlib.h>
#include <functional>
template <typename Func>
void forall(int begin, int end, Func loop_body)
{
  typename std::remove_reference<Func>::type body{loop_body};
  int distance = end-begin;
  #pragma omp target teams distribute parallel for map(to:body)
  for (int idx = 0; idx < distance; ++idx) {
    body(idx);
  }
}
int main() {
  const int n = 10;
  int a[n], b[n], c[n];
  for(int i=0; i<n; ++i) {
    a[i] = 0; b[i] = c[i] = i;
  }
  #pragma omp target enter data map(to: a, b[:n], c)
  forall(0, n, [&] (int i) {
    a[i] += b[i] + c[i];
  });
  #pragma omp target exit data map(from: a[:n]) map(release:b[:n], c[:n])
 
  for(int i=0; i<n; ++i) 
    if (a[i] != 2*i) printf("FAIL\n");
  printf("PASS\n");
  return 0;
}
