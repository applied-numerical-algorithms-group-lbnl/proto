#include <iostream>
#include <Proto.H>

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>

#include <vector>
#include <memory>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "Proto_WriteBoxData.H"
#include "Proto_Timer.H"

#include "test_timer.H"

#if DIM == 3
#define MAX_SIZE_COPY 256
#else
#define MAX_SIZE_COPY 4096
#endif

using namespace std;
using namespace Proto;

typedef Var<double,1> State;
typedef Var<double,1> V;


PROTO_KERNEL_START
unsigned int test_fusion_bc_initialize_stateF(State& a_U, State& a_U2)
{
    a_U(0) = 1;
    a_U2(0) = 0;
    return 0;
}
PROTO_KERNEL_END(test_fusion_bc_initialize_stateF, test_fusion_bc_initialize_state)

void test_fusion_bc_write_data( BoxData<double, 1>&a_state, int it)
{
    char basename[1024];
    sprintf(basename,"euler.%06d",it);

    const char* varnames[1];
    varnames[0] = "martin";
    double origin[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        origin[ii] = 0.0;
    }
    WriteBoxData(basename,a_state,varnames,origin,1);
};

bool test_fusion_bc_check_answer(double *ptr, unsigned int n)
{
#if DIM == 3
  for(int val = 0; val < 2; val++)
  {
    unsigned int shift = val;
    unsigned int value = val + 2;
    for(int x = shift ; x < n - 1 - shift ; x++)
    {
      for(int y = shift ; y < n - 1 - shift ; y++)
      for(int z = shift ; z < n - 1 - shift ; z++)
      {
   	if((y == shift || y == n - 1 - shift) 
	&& (x == shift || x == n - 1 - shift) 
	&& (z == shift || z == n - 1 - shift)) 
          if(ptr[x+y*n+z*n*n] != value)
          {
	    std::cout << " error ("<<x<<","<<y<<","<<z<<"):" <<ptr[x+y*n+z*n*n]<<"!="<< value << std::endl;
            return false;
          }
      }
    }
  }

  for(int x = 2 ; x < n - 3 ; x++)
  {
    for(int y = 2 ; y < n - 3 ; y++)
    for(int z = 2 ; z < n - 3 ; z++)
    {
      if(ptr[x+y*n+z*n*n] != 4)
      {
	 std::cout << " error ("<<x<<","<<y<<","<<z<<"):" <<ptr[x+y*n+z*n*n]<<"!="<< 4 << std::endl;
        return false;
      }
    }
  }
#else
  for(int val = 0; val < 2; val++)
  {
    unsigned int shift = val;
    unsigned int value = val + 2;
    for(int x = shift ; x < n - 1 - shift ; x++)
    {
      for(int y = shift ; y < n - 1 - shift ; y++)
      {
   	if((y == shift || y == n - 1 - shift) 
	&& (x == shift || x == n - 1 - shift)) 
          if(ptr[x+y*n] != value)
          {
	    std::cout << " error ("<<x<<","<<y<<"):" <<ptr[x+y*n]<<"!="<< value << std::endl;
            return false;
          }
      }
    }
  }

  for(int x = 2 ; x < n - 3 ; x++)
  {
    for(int y = 2 ; y < n - 3 ; y++)
    {
      if(ptr[x+y*n] != 4)
      {
	 std::cout << " error ("<<x<<","<<y<<"):" <<ptr[x+y*n]<<"!="<< 4 << std::endl;
        return false;
      }
    }
  }
#endif
  return true;
}

void test_fusion_bc_print(double* ptr, unsigned int n)
{

  if(n>20) return;

#if DIM == 3
  for(int k = 0 ; k < n ; k++)
  {
  for(int j = 0 ; j < n ; j++)
  {
    for(int i = 0 ; i < n ; i++)
    {
      std::cout << ptr[i+j*n+k*n*n] << " ";    
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << std::endl;
  }
#else
  for(int j = 0 ; j < n ; j++)
  {
    for(int i = 0 ; i < n ; i++)
    {
      std::cout << ptr[i+j*n] << " ";    
    }
    std::cout << std::endl;
  }
#endif
}

bool run_test_fusion_bc(unsigned int size1D)
{
  unsigned int size2D= size1D*size1D;

  Stencil<double> ret2 = ((double)(2))*Shift::Zeros();
  Stencil<double> ret3 = ((double)(3))*Shift::Zeros();
  Stencil<double> ret4 = ((double)(4))*Shift::Zeros();

  Box b = Box::Cube(size1D);

  Box bminus= b.grow(-1);   
  Box bminus2= b.grow(-2);   

  BoxData<double,1> myBoxDatain(b);
  BoxData<double,1> myBoxDataout(b);
  BoxData<double,1> myBoxDataoutfused(b);

  std::vector<Stencil<double>> sten;
  std::vector<Box> bx;
  
  for(int i=0;i<2*DIM;i++) sten.push_back(ret2);
  for(int i=0;i<2*DIM;i++) sten.push_back(ret3);

  forallInPlace(test_fusion_bc_initialize_state, b, myBoxDatain, myBoxDataout);

  Box FxInf = b.face(0,Side::Lo); 
  Box FxSup = b.edge(Point::Basis(1),1); 
  Box FyInf = b.face(1,Side::Lo); 
  Box FySup = b.edge(Point::Basis(0),1); 

  bx.push_back(FxInf);
  bx.push_back(FxSup);
  bx.push_back(FyInf);
  bx.push_back(FySup);

  Box FxInfminus = bminus.face(0,Side::Lo); 
  Box FxSupminus = bminus.edge(Point::Basis(1),1); 
  Box FyInfminus = bminus.face(1,Side::Lo); 
  Box FySupminus = bminus.edge(Point::Basis(0),1); 

  bx.push_back(FxInfminus);
  bx.push_back(FxSupminus);
  bx.push_back(FyInfminus);
  bx.push_back(FySupminus);

#if DIM == 3
  Box FzInf = b.edge(Point::Basis(1),2); 
  Box FzSup = b.edge(Point::Basis(0),2); 

  bx.push_back(FzInf);
  bx.push_back(FzSup);

  Box FzInfminus = bminus.edge(Point::Basis(1),2); 
  Box FzSupminus = bminus.edge(Point::Basis(0),2); 

  bx.push_back(FzInfminus);
  bx.push_back(FzSupminus);
#endif

  for( int it = 0 ; it < sten.size() ; it++)
    sten[it].apply(myBoxDatain, myBoxDataout, bx[it], true, 1.0);
  
  ret4.apply(myBoxDatain, myBoxDataout, bminus2, true, 1.0);

  forallInPlace(test_fusion_bc_initialize_state, b, myBoxDatain, myBoxDataoutfused);
  FusedStencil<double> myTest;
  
  //auto& ptr = myTest.getStencils();
  unsigned int nDim = 2; //only 2D// DIM;
  std::vector<Stencil<double>> vec_stencils[nDim];
  std::vector<Box>                vec_boxes[nDim];

  for(int i=0;i<nDim;i++)
    for(int j=0;j<2*DIM;j++)
    {
      vec_boxes[i].push_back(bx[i*4+j]);
    }

  for(int i=0;i<2*DIM;i++) vec_stencils[0].push_back(ret2);
  for(int i=0;i<2*DIM;i++) vec_stencils[1].push_back(ret3);
 
  myTest.define(vec_stencils, vec_boxes, nDim); //copyInfo();

  myTest.cudaApplyFused(myBoxDatain, myBoxDataoutfused, true, 1.0); 
  ret4.apply(myBoxDatain, myBoxDataoutfused, bminus2, true, 1.0);
//  test_fusion_bc_write_boxData(myBoxDataoutfused,2);

  // checking part
#if DIM == 3
  double *h_ptr = new double[size2D*size1D];
  unsigned int nBytes = size2D * size1D * sizeof(double);
#else
  double *h_ptr = new double[size2D];
  unsigned int nBytes = size2D * sizeof(double);
#endif

  double *d_ptr = myBoxDataout.data();
  protoMemcpyGPU(h_ptr,d_ptr,nBytes,protoMemcpyDeviceToHost);
  bool check1 = test_fusion_bc_check_answer(h_ptr,size1D);
//  if(!check1) test_fusion_bc_print(h_ptr,size1D);
  assert(check1);

  d_ptr = myBoxDataoutfused.data();
  protoMemcpyGPU(h_ptr,d_ptr,nBytes,protoMemcpyDeviceToHost);
  bool check2 = test_fusion_bc_check_answer(h_ptr,size1D);
  if(!check2) test_fusion_bc_print(h_ptr,size1D);
  assert(check2);
  return (check1 && check2);
}

bool run_test_fusion_bc_base()
{
  unsigned int size = 8;
  return run_test_fusion_bc(size);
}

bool run_test_fusion_bc_debug()
{
  unsigned int size = MAX_SIZE_COPY;
  return run_test_fusion_bc(size);
}

bool run_test_fusion_bc_stress()
{
        test_timer timer;
	unsigned int size = 8;
	const unsigned int max_size = MAX_SIZE_COPY;
	const unsigned int factor = 2;

        bool b =true;
	while(size<=max_size)
        {
                timer.begin();
                b = run_test_fusion_bc(size);
                timer.end();
                if(b) std::cout << " run_test_irreg_copy_stress_" << size << " ok ... " << timer.duration() << " ms" << std::endl;
                else  std::cout << " error for size = " << size << std::endl;
		if(!b) return false;
                size *= factor;
        }

  return true;
}

bool run_test_fusion_bc_stress_repeated()
{
        test_timer timer;
	unsigned int size = MAX_SIZE_COPY;

        bool b =true;
	for(int i=0 ; i<100 ; i++)
        {
                timer.begin();
                b = run_test_fusion_bc(size);
                timer.end();
                if(b) std::cout << " run_test_irreg_copy_stress_repeated_" << i << " ok ... " << timer.duration() << " ms" << std::endl;
                else  std::cout << " error for size = " << size << std::endl;
		if(!b) return false;
        }

  return true;
}
