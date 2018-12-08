/*
  Copyright Marcin Krotkiewski, University of Oslo, 2012
*/

#define stencil_3x3(c0, c1, c2, shm, tx, ty, bx)			\
  c0*((shm)[tx+0+(ty+0)*bx]) +						\
  c1*((shm)[tx-1+(ty+0)*bx]  + (shm)[tx+1+(ty+0)*bx] + (shm)[tx+0+(ty-1)*bx] + (shm)[tx+0+(ty+1)*bx]) + \
  c2*((shm)[tx-1+(ty-1)*bx]  + (shm)[tx+1+(ty-1)*bx] + (shm)[tx-1+(ty+1)*bx] + (shm)[tx+1+(ty+1)*bx])
  

#define stencil_3x3_reg(c0, c1, c2)					\
  c0*r5 +								\
  c1*(r2+r4+r6+r8) +							\
  c2*(r1+r3+r7+r9)

#define push_regs_exp(shm, bx)			\
  {						\
    r1=(shm)[tx-1+(ty-1)*bx];			\
    r2=(shm)[tx+0+(ty-1)*bx];			\
    r3=(shm)[tx+1+(ty-1)*bx];			\
						\
    r4=(shm)[tx-1+(ty+0)*bx];			\
    r5=(shm)[tx+0+(ty+0)*bx];			\
    r6=(shm)[tx+1+(ty+0)*bx];			\
	    					\
    r7=(shm)[tx-1+(ty+1)*bx];			\
    r8=(shm)[tx+0+(ty+1)*bx];			\
    r9=(shm)[tx+1+(ty+1)*bx];			\
  }						\

__global__ void stencil27_symm_exp_tex(mfloat *out, mfloat a, mfloat b,
				       uint dimx, uint dimy, uint dimz, uint pitch, uint pitchy, uint texoffset, 
				       uint kstart, uint kend)
{
  const uint tx = threadIdx.x;
  const uint ty = threadIdx.y;
  const  int ix = blockIdx.x*blockDim.x + threadIdx.x;	
  const  int iy = blockIdx.y*blockDim.y + threadIdx.y;

  const uint ti = threadIdx.y*blockDim.x + threadIdx.x;
  const uint bxe= blockDim.x+2*8;
  const uint txe= ti%bxe;
  const uint tye= ti/bxe;
  const uint txe2= (ti+blockDim.x*blockDim.y)%bxe;
  const uint tye2= (ti+blockDim.x*blockDim.y)/bxe;
  int  ixe= blockIdx.x*blockDim.x + txe - 8;
  int  iye= blockIdx.y*blockDim.y + tye - 1;
  int  ixe2= blockIdx.x*blockDim.x + txe2 - 8;
  int  iye2= blockIdx.y*blockDim.y + tye2 - 1;
#ifndef MSINGLE
  int2 v;
#endif

  // periodicity
  if(ixe<0)       ixe  += dimx;
  if(ixe>dimx-1)  ixe  -= dimx;
  if(ixe2<0)      ixe2 += dimx;
  if(ixe2>dimx-1) ixe2 -= dimx;

  if(iye<0)       iye  += dimy;
  if(iye>dimy-1)  iye  -= dimy;
  if(iye2<0)      iye2 += dimy;
  if(iye2>dimy-1) iye2 -= dimy;

  mfloat t1 = 0;
  mfloat t2 = 0;
  mfloat t3 = 0;
  mfloat *kernel = d_kernel_3c;
  mfloat C0, C1, C2, C3;
  C0 = kernel[9+4];
  C1 = kernel[4];
  C2 = kernel[1];
  C3 = kernel[0];
  uint i1, i2;

  uint kk;						
  extern __shared__ mfloat shm[];			
  const uint bx = blockDim.x+2*8;				

  i1 = ixe+iye*pitch +texoffset;
  i2 = ixe2+iye2*pitch +texoffset;

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
  shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

  __syncthreads();
  //t1 = convolution_3x3(kernel, shm, tx+8, ty+1, bx);
  t1 = stencil_3x3(C1, C2, C3, shm, tx+8, ty+1, bx);
  __syncthreads();

  i1 += pitch*pitchy;
  i2 += pitch*pitchy;

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
  shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

  __syncthreads();
  //t2 = convolution_3x3(kernel, shm, tx+8, ty+1, bx);
  //t1+= convolution_3x3(kernel+9, shm, tx+8, ty+1, bx);
  t2 = stencil_3x3(C1, C2, C3, shm, tx+8, ty+1, bx);
  t1+= stencil_3x3(C0, C1, C2, shm, tx+8, ty+1, bx);
  __syncthreads();

  for(kk=kstart; kk<kend; kk++){

    __syncthreads();

    i1 += pitch*pitchy;
    i2 += pitch*pitchy;

#ifndef MSINGLE
    v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
    v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
    shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
    shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

    __syncthreads();
    //t3 = convolution_3x3(kernel+18, shm, tx+8, ty+1, bx);
    t3 = stencil_3x3(C1, C2, C3, shm, tx+8, ty+1, bx);

    out[ix + iy*pitch + kk*pitch*pitchy] = t1 + t3;
    t1 = t2 + stencil_3x3(C0, C1, C2, shm, tx+8, ty+1, bx);
    //t1 = t2 + convolution_3x3(kernel+9, shm, tx+8, ty+1, bx);
    t2 = t3;
  }
}


__global__ void stencil27_symm_exp_tex_prefetch(mfloat *out, mfloat a, mfloat b,
						uint dimx, uint dimy, uint dimz, uint pitch, uint pitchy, uint texoffset, 
						uint kstart, uint kend)
{
  mfloat r1, r2, r3, r4, r5, r6, r7, r8, r9;

  const uint tx = threadIdx.x;
  const uint ty = threadIdx.y;
  const  int ix = blockIdx.x*blockDim.x + threadIdx.x;	
  const  int iy = blockIdx.y*blockDim.y + threadIdx.y;

  const uint ti = threadIdx.y*blockDim.x + threadIdx.x;
  const uint bxe= blockDim.x+2*8;
  const uint txe= ti%bxe;
  const uint tye= ti/bxe;
  const uint txe2= (ti+blockDim.x*blockDim.y)%bxe;
  const uint tye2= (ti+blockDim.x*blockDim.y)/bxe;
  int  ixe= blockIdx.x*blockDim.x + txe - 8;
  int  iye= blockIdx.y*blockDim.y + tye - 1;
  int  ixe2= blockIdx.x*blockDim.x + txe2 - 8;
  int  iye2= blockIdx.y*blockDim.y + tye2 - 1;
#ifndef MSINGLE
  int2 v;
#endif

  // periodicity
  if(ixe<0)       ixe  += dimx;
  if(ixe>dimx-1)  ixe  -= dimx;
  if(ixe2<0)      ixe2 += dimx;
  if(ixe2>dimx-1) ixe2 -= dimx;

  if(iye<0)       iye  += dimy;
  if(iye>dimy-1)  iye  -= dimy;
  if(iye2<0)      iye2 += dimy;
  if(iye2>dimy-1) iye2 -= dimy;

  mfloat t1 = 0;
  mfloat t2 = 0;
  mfloat t3 = 0;
  mfloat *kernel = d_kernel_3c;
  mfloat C0, C1, C2, C3;
  C0 = kernel[9+4];
  C1 = kernel[4];
  C2 = kernel[1];
  C3 = kernel[0];
  __syncthreads();

  uint i1, i2;

  uint kk;						
  extern __shared__ mfloat shm[];
  const uint bx = blockDim.x+2*8;				

  i1 = ixe+iye*pitch +texoffset;
  i2 = ixe2+iye2*pitch +texoffset;
#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
  shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

  __syncthreads();  
  push_regs_exp(shm+8+bx, bx);  
  __syncthreads();

  i1 += pitch*pitchy;
  i2 += pitch*pitchy;
#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
  shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

  //t1 = convolution_3x3_reg((kernel));
  t1 = stencil_3x3_reg(C1, C2, C3);

  __syncthreads();  
  push_regs_exp(shm+8+bx, bx);  
  __syncthreads();

  i1 += pitch*pitchy;
  i2 += pitch*pitchy;

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
  shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

  //t2 = convolution_3x3_reg((kernel));
  //t1+= convolution_3x3_reg((kernel+9));
  t2 = stencil_3x3_reg(C1, C2, C3);
  t1+= stencil_3x3_reg(C0, C1, C2);

  for(kk=kstart; kk<kend-1; kk++){

    __syncthreads();  
    push_regs_exp(shm+8+bx, bx);  
    __syncthreads();

    i1 += pitch*pitchy;
    i2 += pitch*pitchy;

#ifndef MSINGLE
    v = tex1Dfetch(texData1D, i1); shm[txe +tye *bx] = __hiloint2double(v.y, v.x);
    v = tex1Dfetch(texData1D, i2); shm[txe2+tye2*bx] = __hiloint2double(v.y, v.x);
#else
    shm[txe +tye *bx] = tex1Dfetch(texData1D, i1);
    shm[txe2+tye2*bx] = tex1Dfetch(texData1D, i2);
#endif

    //t3 = convolution_3x3_reg((kernel+18));
    t3 = stencil_3x3_reg(C1, C2, C3);

    out[ix + iy*pitch + kk*pitch*pitchy] = t1 + t3;
    //t1 = t2 + convolution_3x3_reg((kernel+9));
    t1 = t2 + stencil_3x3_reg(C0, C1, C2);
    t2 = t3;

  }

  __syncthreads();  
  push_regs_exp(shm+8+bx, bx);  
  __syncthreads();

  //*out = t1 + convolution_3x3_reg((kernel+18));
  out[ix + iy*pitch + kk*pitch*pitchy] = t1 + stencil_3x3_reg(C1, C2, C3);
}



__global__ void stencil27_symm_exp_tex_new(mfloat *out, mfloat a, mfloat b,
					   uint dimx, uint dimy, uint dimz, uint pitch, uint pitchy, uint texoffset, 
					   uint kstart, uint kend)
{
  const uint tx = threadIdx.x;
  const uint ty = threadIdx.y;
  const uint ix = blockIdx.x*32 + threadIdx.x;
  const uint iy = blockIdx.y*6  + threadIdx.y;
  const uint ti = threadIdx.y*32 + threadIdx.x;
  const uint tye= ti/48;
  const uint txe= ti-tye*48;
  const uint tye2=tye+4;

  int  ixe = blockIdx.x*32 + txe  - 8;
  int  iye = blockIdx.y*6  + tye  - 1;
  int  iye2= blockIdx.y*6  + tye2 - 1;
#ifndef MSINGLE
  int2 v;
#endif

  // periodicity
  if(ixe<0)       ixe  += dimx;
  if(ixe>dimx-1)  ixe  -= dimx;
  if(iye<0)       iye  += dimy;
  if(iye>dimy-1)  iye  -= dimy;
  if(iye2<0)      iye2 += dimy;
  if(iye2>dimy-1) iye2 -= dimy;

  uint i1, i2;
  
  i1 = ixe+iye*pitch + texoffset;
  i2 = ixe+iye2*pitch + texoffset;

  mfloat t1 = 0;
  mfloat t2 = 0;
  mfloat t3 = 0;
  mfloat *kernel = d_kernel_3c;
  mfloat C0, C1, C2, C3;
  C0 = kernel[9+4];
  C1 = kernel[4];
  C2 = kernel[1];
  C3 = kernel[0];

  uint kk;						
  extern __shared__ mfloat shm[];			

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *48] = tex1Dfetch(texData1D, i1);
  shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

  __syncthreads();
  //t1 = convolution_3x3(kernel, shm, tx+8, ty+1, 48);
  t1 = stencil_3x3(C1, C2, C3, shm, tx+8, ty+1, 48);
  __syncthreads();

  i1 += pitch*pitchy;
  i2 += pitch*pitchy;

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *48] = tex1Dfetch(texData1D, i1);
  shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

  __syncthreads();
  //t2 = convolution_3x3(kernel, shm, tx+8, ty+1, 48);
  //t1+= convolution_3x3(kernel+9, shm, tx+8, ty+1, 48);
  t2 = stencil_3x3(C1, C2, C3, shm, tx+8, ty+1, 48);
  t1+= stencil_3x3(C0, C1, C2, shm, tx+8, ty+1, 48);
  __syncthreads();

  for(kk=kstart; kk<kend; kk++){

    __syncthreads();

    i1 += pitch*pitchy;
    i2 += pitch*pitchy;

#ifndef MSINGLE
    v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
    v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
    shm[txe +tye *48] = tex1Dfetch(texData1D, i1);
    shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

    __syncthreads();
    //t3 = convolution_3x3(kernel+18, shm, tx+8, ty+1, 48);
    t3 = stencil_3x3(C1, C2, C3, shm, tx+8, ty+1, 48);

    out[ix + iy*pitch + kk*pitch*pitchy] = t1 + t3;
    t1 = t2 + stencil_3x3(C0, C1, C2, shm, tx+8, ty+1, 48);
    //t1 = t2 + convolution_3x3(kernel+9, shm, tx+8, ty+1, 48);
    t2 = t3;
  }
}


__global__ void stencil27_symm_exp_tex_prefetch_new(mfloat *out, mfloat a, mfloat b,
						    uint dimx, uint dimy, uint dimz, uint pitch, uint pitchy, uint texoffset, 
						    uint kstart, uint kend)
{
  mfloat r1, r2, r3, r4, r5, r6, r7, r8, r9;

  const uint tx = threadIdx.x;
  const uint ty = threadIdx.y;
  const uint ix = blockIdx.x*32 + threadIdx.x;
  const uint iy = blockIdx.y*6  + threadIdx.y;
  const uint ti = threadIdx.y*32 + threadIdx.x;
  const uint tye= ti/48;
  const uint txe= ti-tye*48;
  const uint tye2=tye+4;

  int  ixe = blockIdx.x*32 + txe  - 8;
  int  iye = blockIdx.y*6  + tye  - 1;
  int  iye2= blockIdx.y*6  + tye2 - 1;
#ifndef MSINGLE
  int2 v;
#endif

  // periodicity
  if(ixe<0)       ixe  += dimx;
  if(ixe>dimx-1)  ixe  -= dimx;
  if(iye<0)       iye  += dimy;
  if(iye>dimy-1)  iye  -= dimy;
  if(iye2<0)      iye2 += dimy;
  if(iye2>dimy-1) iye2 -= dimy;

  uint i1, i2;
  
  i1 = ixe+iye*pitch + texoffset;
  i2 = ixe+iye2*pitch + texoffset;

  mfloat t1 = 0;
  mfloat t2 = 0;
  mfloat t3 = 0;
  mfloat *kernel = d_kernel_3c;
  mfloat C0, C1, C2, C3;
  C0 = kernel[9+4];
  C1 = kernel[4];
  C2 = kernel[1];
  C3 = kernel[0];

  __syncthreads();

  uint kk;						
  extern __shared__ mfloat shm[];

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *48] = tex1Dfetch(texData1D, i1);
  shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

  __syncthreads();  
  push_regs_exp(shm+8+48, 48);  
  __syncthreads();

  i1 += pitch*pitchy;
  i2 += pitch*pitchy;

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye *48] = tex1Dfetch(texData1D, i1);
  shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

  //t1 = convolution_3x3_reg((kernel));
  t1 = stencil_3x3_reg(C1, C2, C3);

  __syncthreads();  
  push_regs_exp(shm+8+48, 48);  
  __syncthreads();

  i1 += pitch*pitchy;
  i2 += pitch*pitchy;

#ifndef MSINGLE
  v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
  v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
  shm[txe +tye*48] = tex1Dfetch(texData1D, i1);
  shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

  //t2 = convolution_3x3_reg((kernel));
  //t1+= convolution_3x3_reg((kernel+9));
  t2 = stencil_3x3_reg(C1, C2, C3);
  t1+= stencil_3x3_reg(C0, C1, C2);

  for(kk=kstart; kk<kend-1; kk++){ // 1

    __syncthreads();  
    push_regs_exp(shm+8+48, 48);  
    __syncthreads();

    i1 += pitch*pitchy;
    i2 += pitch*pitchy;

#ifndef MSINGLE
    v = tex1Dfetch(texData1D, i1); shm[txe +tye *48] = __hiloint2double(v.y, v.x);
    v = tex1Dfetch(texData1D, i2); shm[txe+tye2*48] = __hiloint2double(v.y, v.x);
#else
    shm[txe +tye *48] = tex1Dfetch(texData1D, i1);
    shm[txe+tye2*48] = tex1Dfetch(texData1D, i2);
#endif

    //t3 = convolution_3x3_reg((kernel+18));
    // 13
    t3 = stencil_3x3_reg(C1, C2, C3);

    out[ix + iy*pitch + kk*pitch*pitchy] = t1 + t3; // 2
    //t1 = t2 + convolution_3x3_reg((kernel+9));
    t1 = t2 + stencil_3x3_reg(C0, C1, C2);
    t2 = t3; // 3
  }

  __syncthreads();  
  push_regs_exp(shm+8+48, 48);  
  __syncthreads();

  //*out = t1 + convolution_3x3_reg((kernel+18));
  out[ix + iy*pitch + kk*pitch*pitchy] = t1 + stencil_3x3_reg(C1, C2, C3);
}
