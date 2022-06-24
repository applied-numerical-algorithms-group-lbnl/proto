#include "Proto.H"

int main()
{
#ifdef PROTO_CUDA
  protoSetDevice(5);
#endif
  unsigned int size = 8;
  unsigned int nbAllocation = 10;
  double ** device = new double*[nbAllocation];
  double ** host = new double*[nbAllocation];

  unsigned int acc = 0;

  for(int i = 0 ; i < nbAllocation ; i++)
  {
    double *tmp;
    protoMalloc(Proto::MEMTYPE_DEFAULT, tmp, size*sizeof(double)*i);
    device[i] = tmp;
    acc += size*sizeof(double)*i;
  }

  std::cout << " first cumulate malloc size: " << acc << std::endl;

  for(int i = 0 ; i < nbAllocation ; i++)
    protoFree(Proto::MEMTYPE_DEFAULT,device[i]);

  acc = 0;

  for(int i = 0 ; i < nbAllocation ; i+=3)
  {
    double *tmp;
    protoMalloc(Proto::HOST, tmp, size*sizeof(double)*i);
    host[i] = tmp;
    acc += size*sizeof(double)*i;
  }

  std::cout << " first cumulate malloc size on host: " << acc << std::endl;

  for(int i = 0 ; i < nbAllocation ; i+=2)
    protoFree(Proto::HOST,host[i]);

  acc = 0;

  for(int i = 0 ; i < nbAllocation ; i+=2)
  {
    double *tmp;
    protoMalloc(Proto::MEMTYPE_DEFAULT, tmp, size*sizeof(double)*i);
    device[i] = tmp;
    acc += size*sizeof(double)*i;
  }

  std::cout << " second cumulate malloc size: " << acc << std::endl;

  for(int i = 0 ; i < nbAllocation ; i+=2)
    protoFree(Proto::MEMTYPE_DEFAULT,device[i]);

  // PRINT_MEMORY_INFO(); macro toggled off in Proto_MemInfo.H
}
