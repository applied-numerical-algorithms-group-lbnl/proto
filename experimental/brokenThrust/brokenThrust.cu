#include <thrust/device_vector.h>

namespace Bob
{

  thrust::device_vector<double> testObject;

}

int main(int argc, char* argv[])
{

  Bob::testObject.resize(25);

  Bob::testObject.clear();
  return 0;
}
