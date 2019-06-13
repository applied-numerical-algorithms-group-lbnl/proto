#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define min(x,y) (((x)<(y))?(x):(y))
#define max(x,y) (((x)>(y))?(x):(y))
#define abs(x) ((x)<0?-(x):(x))
#define absmin(x,y) ((x)=min(abs((x)),abs((y))))
#define absmax(x,y) ((x)=max(abs((x)),abs((y))))
#define floord(x,y) ((x)/(y))
#define sgn(x) ((x)<0?-1:1)
#define offset2(i,j,M) ((j)+(i)*(M))
#define offset3(i,j,k,M,N) ((k)+((j)+(i)*(M))*(N))
#define offset4(i,j,k,l,M,N,P) ((l)+((k)+((j)+(i)*(M))*(N))*(P))
#define arrinit(ptr,val,size) for(unsigned __i__=0;__i__<(size);__i__++) (ptr)[__i__]=(val)
#define arrprnt(name,arr,size) {\
fprintf(stderr,"%s={",(name));\
for(unsigned __i__=0;__i__<(size);__i__++) fprintf(stderr,"%lg,",(arr)[__i__]);\
fprintf(stderr,"}\n");}
#define F_ave_f_d1(c,y,x) F_ave_f_d1[offset3((c),(y)+2,(x),((17+2+1)),(16+1))]
#define F_ave_f_d2(c,y,x) F_ave_f_d2[offset3((c),(y),(x)+2,((16+1)),(17+2+1))]
#define F_bar_f_d1(c,y,x) F_bar_f_d1[offset3((c),(y)+3,(x),((18+3+1)),(16+1))]
#define F_bar_f_d2(c,y,x) F_bar_f_d2[offset3((c),(y),(x)+3,((16+1)),(18+3+1))]
#define F_div_f_d1(c,y,x) F_div_f_d1[offset3((c),(y)+2,(x),((17+2+1)),(15+1))]
#define F_div_f_d2(c,y,x) F_div_f_d2[offset3((c),(y),(x)+2,((15+1)),(17+2+1))]
#define F_lap_f_d1(c,y,x) F_lap_f_d1[offset3((c),(y)+2,(x),((17+2+1)),(16+1))]
#define F_lap_f_d2(c,y,x) F_lap_f_d2[offset3((c),(y),(x)+2,((16+1)),(17+2+1))]
#define U(c,y,x) U[offset3((c),(y)+4,(x)+4,((19+4+1)),(19+4+1))]
#define W(c,y,x) W[offset3((c),(y)+3,(x)+3,((18+3+1)),(18+3+1))]
#define W_ave(c,y,x) W_ave[offset3((c),(y)+3,(x)+3,((18+3+1)),(18+3+1))]
#define W_aveH_d1(c,y,x) W_aveH_d1[offset3((c),(y)+3,(x)+1,((18+3+1)),(16+1+1))]
#define W_aveH_d2(c,y,x) W_aveH_d2[offset3((c),(y)+1,(x)+3,((16+1+1)),(18+3+1))]
#define W_aveL_d1(c,y,x) W_aveL_d1[offset3((c),(y)+3,(x),((18+3+1)),(17+1))]
#define W_aveL_d2(c,y,x) W_aveL_d2[offset3((c),(y),(x)+3,((17+1)),(18+3+1))]
#define W_ave_f_d1(c,y,x) W_ave_f_d1[offset3((c),(y)+3,(x),((18+3+1)),(16+1))]
#define W_ave_f_d2(c,y,x) W_ave_f_d2[offset3((c),(y),(x)+3,((16+1)),(18+3+1))]
#define W_bar(c,y,x) W_bar[offset3((c),(y)+4,(x)+4,((19+4+1)),(19+4+1))]
#define W_f_d1(c,y,x) W_f_d1[offset3((c),(y)+2,(x),((17+2+1)),(16+1))]
#define W_f_d2(c,y,x) W_f_d2[offset3((c),(y),(x)+2,((16+1)),(17+2+1))]
#define rhs(c,y,x) rhs[offset3((c),(y),(x),((15+1)),(15+1))]
#define u(c,y,x) u[offset3((c),(y)+3,(x)+3,((18+3+1)),(18+3+1))]
#define umax(y,x) umax[offset2((y),(x),(15+1))]

double euler_step(const double* U, double* rhs);
inline double euler_step(const double* U, double* rhs) {
    int t1,t2,t3,t4,t5;
    double* W_bar = (double*) calloc(((4)*(19+4+1))*(19+4+1),sizeof(double));
    double* u = (double*) calloc(((4)*(18+3+1))*(18+3+1),sizeof(double));
    double* W = (double*) calloc(((4)*(18+3+1))*(18+3+1),sizeof(double));
    double* umax = (double*) calloc(((15+1))*(15+1),sizeof(double));
    double retval = 0;
    double* W_ave = (double*) calloc(((4)*(18+3+1))*(18+3+1),sizeof(double));
    double* W_aveL_d1 = (double*) calloc(((4)*(18+3+1))*(17+1),sizeof(double));
    double* W_aveH_d1 = (double*) calloc(((4)*(18+3+1))*(16+1+1),sizeof(double));
    double* W_ave_f_d1 = (double*) calloc(((4)*(18+3+1))*(16+1),sizeof(double));
    double* F_bar_f_d1 = (double*) calloc(((4)*(18+3+1))*(16+1),sizeof(double));
    double* W_f_d1 = (double*) calloc(((4)*(17+2+1))*(16+1),sizeof(double));
    double* F_ave_f_d1 = (double*) calloc(((4)*(17+2+1))*(16+1),sizeof(double));
    double* F_lap_f_d1 = (double*) calloc(((4)*(17+2+1))*(16+1),sizeof(double));
    double* F_div_f_d1 = (double*) calloc(((4)*(17+2+1))*(15+1),sizeof(double));
    double* W_aveL_d2 = (double*) calloc(((4)*(17+1))*(18+3+1),sizeof(double));
    double* W_aveH_d2 = (double*) calloc(((4)*(16+1+1))*(18+3+1),sizeof(double));
    double* W_ave_f_d2 = (double*) calloc(((4)*(16+1))*(18+3+1),sizeof(double));
    double* F_bar_f_d2 = (double*) calloc(((4)*(16+1))*(18+3+1),sizeof(double));
    double* W_f_d2 = (double*) calloc(((4)*(16+1))*(17+2+1),sizeof(double));
    double* F_ave_f_d2 = (double*) calloc(((4)*(16+1))*(17+2+1),sizeof(double));
    double* F_lap_f_d2 = (double*) calloc(((4)*(16+1))*(17+2+1),sizeof(double));
    double* F_div_f_d2 = (double*) calloc(((4)*(15+1))*(17+2+1),sizeof(double));

// consToPrim1
#undef s0
#define s0(y,x) {\
W_bar(0,(y),(x)) = U(0,(y),(x));\
W_bar(1,(y),(x)) = U(1,(y),(x))/U(0,(y),(x));\
W_bar(2,(y),(x)) = U(2,(y),(x))/U(0,(y),(x));\
W_bar(3,(y),(x)) = (U(3,(y),(x))-0.500000*U(0,(y),(x))*(((U(1,(y),(x))/U(0,(y),(x)))*(U(1,(y),(x))/U(0,(y),(x))))+((U(2,(y),(x))/U(0,(y),(x)))*(U(2,(y),(x))/U(0,(y),(x))))))*(1.400000-1.000000);\
}

for(t1 = -4; t1 <= 19; t1++) {
  for(t2 = -4; t2 <= 19; t2++) {
    s0(t1,t2);
  }
}

// deconvolve
#undef s0
#define s0(c,y,x) u((c),(y),(x))=(1.166667*U((c),(y),(x)))+((-0.041667)*U((c),(y),(x)+1))+((-0.041667)*U((c),(y),(x)-1))+((-0.041667)*U((c),(y)+1,(x)))+((-0.041667)*U((c),(y)-1,(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    for(t3 = -3; t3 <= 18; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// consToPrim2
#undef s0
#define s0(y,x) {\
W(0,(y),(x)) = u(0,(y),(x));\
W(1,(y),(x)) = u(1,(y),(x))/u(0,(y),(x));\
W(2,(y),(x)) = u(2,(y),(x))/u(0,(y),(x));\
W(3,(y),(x)) = (u(3,(y),(x))-0.500000*u(0,(y),(x))*(((u(1,(y),(x))/u(0,(y),(x)))*(u(1,(y),(x))/u(0,(y),(x))))+((u(2,(y),(x))/u(0,(y),(x)))*(u(2,(y),(x))/u(0,(y),(x))))))*(1.400000-1.000000);\
}

for(t1 = -3; t1 <= 18; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    s0(t1,t2);
  }
}

// waveSpeedBound1
#undef s0
#define s0(y,x) {\
umax((y),(x))=(2.000000*sqrt(1.400000*W(3,(y),(x))/W(0,(y),(x))))+W(1,(y),(x))+W(2,(y),(x));\
}

for(t1 = 0; t1 <= 15; t1++) {
  for(t2 = 0; t2 <= 15; t2++) {
    s0(t1,t2);
  }
}

// absMax
#undef s0
#define s0(y,x) absmax(retval,umax((y),(x)))

for(t1 = 0; t1 <= 15; t1++) {
  for(t2 = 0; t2 <= 15; t2++) {
    s0(t1,t2);
  }
}

// laplacian
#undef s0
#define s0(c,y,x) W_ave((c),(y),(x))=((-0.166667)*W_bar((c),(y),(x)))+(0.041667*W_bar((c),(y),(x)+1))+(0.041667*W_bar((c),(y),(x)-1))+(0.041667*W_bar((c),(y)+1,(x)))+(0.041667*W_bar((c),(y)-1,(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    for(t3 = -3; t3 <= 18; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// increment
#undef s0
#define s0(c,y,x) W_ave((c),(y),(x))+=W((c),(y),(x))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    for(t3 = -3; t3 <= 18; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// interpL_d1
#undef s0
#define s0(c,y,x) W_aveL_d1((c),(y),(x))=(0.033333*W_ave((c),(y),(x)-3))+((-0.050000)*W_ave((c),(y),(x)+1))+((-0.216667)*W_ave((c),(y),(x)-2))+(0.450000*W_ave((c),(y),(x)))+(0.783333*W_ave((c),(y),(x)-1))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    for(t3 = 0; t3 <= 17; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// interpH_d1
#undef s0
#define s0(c,y,x) W_aveH_d1((c),(y),(x))=((-0.050000)*W_ave((c),(y),(x)-2))+(0.033333*W_ave((c),(y),(x)+2))+(0.450000*W_ave((c),(y),(x)-1))+((-0.216667)*W_ave((c),(y),(x)+1))+(0.783333*W_ave((c),(y),(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    for(t3 = -1; t3 <= 16; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// upwindState1
#undef s0
#define s0(y,x) {\
W_ave_f_d1(0,(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))) ? (W_aveL_d1(0,(y),(x))) : (W_aveH_d1(0,(y),(x)));\
W_ave_f_d1(1,(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))) ? (W_aveL_d1(1,(y),(x))) : (W_aveH_d1(1,(y),(x)));\
W_ave_f_d1(2,(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))) ? (W_aveL_d1(2,(y),(x))) : (W_aveH_d1(2,(y),(x)));\
W_ave_f_d1(3,(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))) ? (W_aveL_d1(3,(y),(x))) : (W_aveH_d1(3,(y),(x)));\
W_ave_f_d1(0,(y),(x))+=(0.000000 < sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))+(-sgn(((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))))*((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000) ? (((((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000+((((((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))*(W_aveL_d1(0+1,(y),(x))-W_aveH_d1(0+1,(y),(x))))*0.500000-W_ave_f_d1(3,(y),(x))))/((sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000)))*sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000)))) : (0.000000);\
W_ave_f_d1(0+1,(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))+(-sgn(((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))))*((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000) ? (((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))) : (W_ave_f_d1(0+1,(y),(x)));\
W_ave_f_d1(3,(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))+(-sgn(((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000+((W_aveL_d1(3,(y),(x))-W_aveH_d1(3,(y),(x))))/((2.000000*(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))))*((W_aveL_d1(0+1,(y),(x))+W_aveH_d1(0+1,(y),(x))))*0.500000) ? (((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000+((((((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(3,(y),(x))+W_aveH_d1(3,(y),(x))))*0.500000))/(((W_aveL_d1(0,(y),(x))+W_aveH_d1(0,(y),(x))))*0.500000))))*(W_aveL_d1(0+1,(y),(x))-W_aveH_d1(0+1,(y),(x))))*0.500000) : (W_ave_f_d1(3,(y),(x)));\
}

for(t1 = -3; t1 <= 18; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    s0(t1,t2);
  }
}

// getFlux1
#undef s0
#define s0(y,x) {\
F_bar_f_d1(0,(y),(x))=(W_ave_f_d1(0+1,(y),(x)))*W_ave_f_d1(0,(y),(x));\
F_bar_f_d1(1,(y),(x))=W_ave_f_d1(1,(y),(x))*F_bar_f_d1(0,(y),(x));\
F_bar_f_d1(2,(y),(x))=W_ave_f_d1(2,(y),(x))*F_bar_f_d1(0,(y),(x));\
F_bar_f_d1(0+1,(y),(x))+=W_ave_f_d1(3,(y),(x));\
F_bar_f_d1(3,(y),(x))=((1.400000/(1.400000-1))*W_ave_f_d1(0+1,(y),(x)))*W_ave_f_d1(3,(y),(x))+0.500000*F_bar_f_d1(0,(y),(x))*((W_ave_f_d1(1,(y),(x))*W_ave_f_d1(1,(y),(x)))+(W_ave_f_d1(2,(y),(x))*W_ave_f_d1(2,(y),(x))));\
}

for(t1 = -3; t1 <= 18; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    s0(t1,t2);
  }
}

// deconvolve_f_d1
#undef s0
#define s0(c,y,x) W_f_d1((c),(y),(x))=(1.083333*W_ave_f_d1((c),(y),(x)))+((-0.041667)*W_ave_f_d1((c),(y)+1,(x)))+((-0.041667)*W_ave_f_d1((c),(y)-1,(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -2; t2 <= 17; t2++) {
    for(t3 = 0; t3 <= 16; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// getFlux2
#undef s0
#define s0(y,x) {\
F_ave_f_d1(0,(y),(x))=(W_f_d1(0+1,(y),(x)))*W_f_d1(0,(y),(x));\
F_ave_f_d1(1,(y),(x))=W_f_d1(1,(y),(x))*F_ave_f_d1(0,(y),(x));\
F_ave_f_d1(2,(y),(x))=W_f_d1(2,(y),(x))*F_ave_f_d1(0,(y),(x));\
F_ave_f_d1(0+1,(y),(x))+=W_f_d1(3,(y),(x));\
F_ave_f_d1(3,(y),(x))=((1.400000/(1.400000-1))*W_f_d1(0+1,(y),(x)))*W_f_d1(3,(y),(x))+0.500000*F_ave_f_d1(0,(y),(x))*((W_f_d1(1,(y),(x))*W_f_d1(1,(y),(x)))+(W_f_d1(2,(y),(x))*W_f_d1(2,(y),(x))));\
}

for(t1 = -2; t1 <= 17; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    s0(t1,t2);
  }
}

// smul_d1
#undef s0
#define s0(c,y,x) F_bar_f_d1((c),(y),(x))*=0.041667

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    for(t3 = 0; t3 <= 16; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// lap_f_d1
#undef s0
#define s0(c,y,x) F_lap_f_d1((c),(y),(x))=((-0.083333)*F_bar_f_d1((c),(y),(x)))+(0.041667*F_bar_f_d1((c),(y)+1,(x)))+(0.041667*F_bar_f_d1((c),(y)-1,(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -2; t2 <= 17; t2++) {
    for(t3 = 0; t3 <= 16; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// inc_f_d1
#undef s0
#define s0(c,y,x) F_ave_f_d1((c),(y),(x))+=F_lap_f_d1((c),(y),(x))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -2; t2 <= 17; t2++) {
    for(t3 = 0; t3 <= 16; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// div_f_d1
#undef s0
#define s0(c,y,x) F_div_f_d1((c),(y),(x))=((-1.000000)*F_ave_f_d1((c),(y),(x)))+(1.000000*F_ave_f_d1((c),(y),(x)+1))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -2; t2 <= 17; t2++) {
    for(t3 = 0; t3 <= 15; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// inc_rhs_d1
#undef s0
#define s0(c,y,x) rhs((c),(y),(x))+=F_div_f_d1((c),(y),(x))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 15; t2++) {
    for(t3 = 0; t3 <= 15; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// interpL_d2
#undef s0
#define s0(c,y,x) W_aveL_d2((c),(y),(x))=(0.033333*W_ave((c),(y)-3,(x)))+((-0.050000)*W_ave((c),(y)+1,(x)))+((-0.216667)*W_ave((c),(y)-2,(x)))+(0.450000*W_ave((c),(y),(x)))+(0.783333*W_ave((c),(y)-1,(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 17; t2++) {
    for(t3 = -3; t3 <= 18; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// interpH_d2
#undef s0
#define s0(c,y,x) W_aveH_d2((c),(y),(x))=((-0.050000)*W_ave((c),(y)-2,(x)))+(0.033333*W_ave((c),(y)+2,(x)))+(0.450000*W_ave((c),(y)-1,(x)))+((-0.216667)*W_ave((c),(y)+1,(x)))+(0.783333*W_ave((c),(y),(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = -1; t2 <= 16; t2++) {
    for(t3 = -3; t3 <= 18; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// upwindState2
#undef s0
#define s0(y,x) {\
W_ave_f_d2(0,(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))) ? (W_aveL_d2(0,(y),(x))) : (W_aveH_d2(0,(y),(x)));\
W_ave_f_d2(1,(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))) ? (W_aveL_d2(1,(y),(x))) : (W_aveH_d2(1,(y),(x)));\
W_ave_f_d2(2,(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))) ? (W_aveL_d2(2,(y),(x))) : (W_aveH_d2(2,(y),(x)));\
W_ave_f_d2(3,(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))) ? (W_aveL_d2(3,(y),(x))) : (W_aveH_d2(3,(y),(x)));\
W_ave_f_d2(0,(y),(x))+=(0.000000 < sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))+(-sgn(((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))))*((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000) ? (((((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000+((((((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))*(W_aveL_d2(1+1,(y),(x))-W_aveH_d2(1+1,(y),(x))))*0.500000-W_ave_f_d2(3,(y),(x))))/((sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000)))*sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000)))) : (0.000000);\
W_ave_f_d2(1+1,(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))+(-sgn(((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))))*((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000) ? (((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))) : (W_ave_f_d2(1+1,(y),(x)));\
W_ave_f_d2(3,(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))+(-sgn(((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000+((W_aveL_d2(3,(y),(x))-W_aveH_d2(3,(y),(x))))/((2.000000*(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))))*((W_aveL_d2(1+1,(y),(x))+W_aveH_d2(1+1,(y),(x))))*0.500000) ? (((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000+((((((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(3,(y),(x))+W_aveH_d2(3,(y),(x))))*0.500000))/(((W_aveL_d2(0,(y),(x))+W_aveH_d2(0,(y),(x))))*0.500000))))*(W_aveL_d2(1+1,(y),(x))-W_aveH_d2(1+1,(y),(x))))*0.500000) : (W_ave_f_d2(3,(y),(x)));\
}

for(t1 = 0; t1 <= 16; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    s0(t1,t2);
  }
}

// getFlux3
#undef s0
#define s0(y,x) {\
F_bar_f_d2(0,(y),(x))=(W_ave_f_d2(1+1,(y),(x)))*W_ave_f_d2(0,(y),(x));\
F_bar_f_d2(1,(y),(x))=W_ave_f_d2(1,(y),(x))*F_bar_f_d2(0,(y),(x));\
F_bar_f_d2(2,(y),(x))=W_ave_f_d2(2,(y),(x))*F_bar_f_d2(0,(y),(x));\
F_bar_f_d2(1+1,(y),(x))+=W_ave_f_d2(3,(y),(x));\
F_bar_f_d2(3,(y),(x))=((1.400000/(1.400000-1))*W_ave_f_d2(1+1,(y),(x)))*W_ave_f_d2(3,(y),(x))+0.500000*F_bar_f_d2(0,(y),(x))*((W_ave_f_d2(1,(y),(x))*W_ave_f_d2(1,(y),(x)))+(W_ave_f_d2(2,(y),(x))*W_ave_f_d2(2,(y),(x))));\
}

for(t1 = 0; t1 <= 16; t1++) {
  for(t2 = -3; t2 <= 18; t2++) {
    s0(t1,t2);
  }
}

// deconvolve_f_d2
#undef s0
#define s0(c,y,x) W_f_d2((c),(y),(x))=(1.083333*W_ave_f_d2((c),(y),(x)))+((-0.041667)*W_ave_f_d2((c),(y),(x)+1))+((-0.041667)*W_ave_f_d2((c),(y),(x)-1))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    for(t3 = -2; t3 <= 17; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// getFlux4
#undef s0
#define s0(y,x) {\
F_ave_f_d2(0,(y),(x))=(W_f_d2(1+1,(y),(x)))*W_f_d2(0,(y),(x));\
F_ave_f_d2(1,(y),(x))=W_f_d2(1,(y),(x))*F_ave_f_d2(0,(y),(x));\
F_ave_f_d2(2,(y),(x))=W_f_d2(2,(y),(x))*F_ave_f_d2(0,(y),(x));\
F_ave_f_d2(1+1,(y),(x))+=W_f_d2(3,(y),(x));\
F_ave_f_d2(3,(y),(x))=((1.400000/(1.400000-1))*W_f_d2(1+1,(y),(x)))*W_f_d2(3,(y),(x))+0.500000*F_ave_f_d2(0,(y),(x))*((W_f_d2(1,(y),(x))*W_f_d2(1,(y),(x)))+(W_f_d2(2,(y),(x))*W_f_d2(2,(y),(x))));\
}

for(t1 = 0; t1 <= 16; t1++) {
  for(t2 = -2; t2 <= 17; t2++) {
    s0(t1,t2);
  }
}

// smul_d2
#undef s0
#define s0(c,y,x) F_bar_f_d2((c),(y),(x))*=0.041667

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    for(t3 = -3; t3 <= 18; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// lap_f_d2
#undef s0
#define s0(c,y,x) F_lap_f_d2((c),(y),(x))=((-0.083333)*F_bar_f_d2((c),(y),(x)))+(0.041667*F_bar_f_d2((c),(y),(x)+1))+(0.041667*F_bar_f_d2((c),(y),(x)-1))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    for(t3 = -2; t3 <= 17; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// inc_f_d2
#undef s0
#define s0(c,y,x) F_ave_f_d2((c),(y),(x))+=F_lap_f_d2((c),(y),(x))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 16; t2++) {
    for(t3 = -2; t3 <= 17; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// div_f_d2
#undef s0
#define s0(c,y,x) F_div_f_d2((c),(y),(x))=((-1.000000)*F_ave_f_d2((c),(y),(x)))+(1.000000*F_ave_f_d2((c),(y)+1,(x)))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 15; t2++) {
    for(t3 = -2; t3 <= 17; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// inc_rhs_d2
#undef s0
#define s0(c,y,x) rhs((c),(y),(x))+=F_div_f_d2((c),(y),(x))

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 15; t2++) {
    for(t3 = 0; t3 <= 15; t3++) {
      s0(t1,t2,t3);
    }
  }
}

// muldx
#undef s0
#define s0(c,y,x) rhs((c),(y),(x))*=-1.000000

for(t1 = 0; t1 <= 3; t1++) {
  for(t2 = 0; t2 <= 15; t2++) {
    for(t3 = 0; t3 <= 15; t3++) {
      s0(t1,t2,t3);
    }
  }
}

    free(W_bar);
    free(u);
    free(W);
    free(umax);
    free(W_ave);
    free(W_aveL_d1);
    free(W_aveH_d1);
    free(W_ave_f_d1);
    free(F_bar_f_d1);
    free(W_f_d1);
    free(F_ave_f_d1);
    free(F_lap_f_d1);
    free(F_div_f_d1);
    free(W_aveL_d2);
    free(W_aveH_d2);
    free(W_ave_f_d2);
    free(F_bar_f_d2);
    free(W_f_d2);
    free(F_ave_f_d2);
    free(F_lap_f_d2);
    free(F_div_f_d2);

    return (retval);
}    // euler_step

#undef min
#undef max
#undef abs
#undef absmin
#undef absmax
#undef floord
#undef sgn
#undef offset2
#undef offset3
#undef offset4
#undef arrinit
#undef arrprnt
#undef F_ave_f_d1
#undef F_ave_f_d2
#undef F_bar_f_d1
#undef F_bar_f_d2
#undef F_div_f_d1
#undef F_div_f_d2
#undef F_lap_f_d1
#undef F_lap_f_d2
#undef U
#undef W
#undef W_ave
#undef W_aveH_d1
#undef W_aveH_d2
#undef W_aveL_d1
#undef W_aveL_d2
#undef W_ave_f_d1
#undef W_ave_f_d2
#undef W_bar
#undef W_f_d1
#undef W_f_d2
#undef rhs
#undef u
#undef umax

