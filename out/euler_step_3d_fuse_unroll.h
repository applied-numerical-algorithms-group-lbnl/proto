// 'euler_step' code generated by 'edavis' at 06/26/2019 11:28:36
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
#define U(c,z,y,x) (((z) >= -4 && (z) <= 67 && (y) >= -4 && (y) <= 67 && (x) >= -4 && (x) <= 67) ? U[offset4((c),max(0,(z)+4),max(0,(y)+4),max(0,(x)+4),(67+4+1),(67+4+1),(67+4+1))] : 0.0)
#define W_bar(c,z,y,x) W_bar[offset4((c),(z)+5,(y)+5,(x)+4,(7),(7),(5)) & 2047]
#define u(c,z,y,x) u[(c)]
#define W(c,z,y,x) W[(c)]
#define umax(z,y,x) umax
#define W_ave(c,z,y,x) W_ave[offset3((z)+7,(y)+7,(x)+7,(10),(10)) & 1023]
#define W_aveL_d1(c,z,y,x) W_aveL_d1[(c)]
#define W_aveH_d1(c,z,y,x) W_aveH_d1[(c)]
#define W_ave_f_d1(c,z,y,x) W_ave_f_d1[offset4((c),(z)+5,(y)+5,(x)+4,(7),(7),(5)) & 2047]
#define F_bar_f_d1(c,z,y,x) F_bar_f_d1[offset4((c),(z)+5,(y)+5,(x)+4,(7),(7),(5)) & 2047]
#define W_f_d1(c,z,y,x) W_f_d1[(c)]
#define F_ave_f_d1(c,z,y,x) F_ave_f_d1[offset2((c),(x)+4,6) & 31]
#define F_lap_f_d1(c,z,y,x) F_lap_f_d1
#define F_div_f_d1(c,z,y,x) F_div_f_d1
#define rhs(c,z,y,x) rhs[offset4((c),(z),(y),(x),(63+1),(63+1),(63+1))]
#define W_aveL_d2(c,z,y,x) W_aveL_d2[(c)]
#define W_aveH_d2(c,z,y,x) W_aveH_d2[(c)]
#define W_ave_f_d2(c,z,y,x) W_ave_f_d2[offset4((c),(z)+5,(y)+4,(x)+5,(7),(5),(7)) & 2047]
#define F_bar_f_d2(c,z,y,x) F_bar_f_d2[offset4((c),(z)+5,(y)+4,(x)+5,(7),(5),(7)) & 2047]
#define W_f_d2(c,z,y,x) W_f_d2[(c)]
#define F_ave_f_d2(c,z,y,x) F_ave_f_d2[offset2((c),(y)+4,6) & 31]
#define F_lap_f_d2(c,z,y,x) F_lap_f_d2
#define F_div_f_d2(c,z,y,x) F_div_f_d2
#define W_aveL_d3(c,z,y,x) W_aveL_d3[(c)]
#define W_aveH_d3(c,z,y,x) W_aveH_d3[(c)]
#define W_ave_f_d3(c,z,y,x) W_ave_f_d3[offset4((c),(z)+4,(y)+5,(x)+5,(5),(7),(7)) & 2047]
#define F_bar_f_d3(c,z,y,x) F_bar_f_d3[offset4((c),(z)+4,(y)+5,(x)+5,(5),(7),(7)) & 2047]
#define W_f_d3(c,z,y,x) W_f_d3[(c)]
#define F_ave_f_d3(c,z,y,x) F_ave_f_d3[offset2((c),(z)+4,6) & 31]
#define F_lap_f_d3(c,z,y,x) F_lap_f_d3
#define F_div_f_d3(c,z,y,x) F_div_f_d3

double euler_step(const double* U, double* rhs);
inline double euler_step(const double* U, double* rhs) {
    int t1,t2,t3,t4,t5;
    int z,y,x;
    double rho, px, py, pz, e;
    double W_bar[2048];
    double u[5];
    double W[5];
    double umax;
    double retval;
    double W_ave[1024]; //10*10*10];
    double W_aveL_d1[5];
    double W_aveH_d1[5];
    double W_ave_f_d1[2048];
    double F_bar_f_d1[2048];
    double W_f_d1[5];
    double F_ave_f_d1[32];
    double F_lap_f_d1;
    double F_div_f_d1;
    double W_aveL_d2[5];
    double W_aveH_d2[5];
    double W_ave_f_d2[2048];
    double F_bar_f_d2[2048];
    double W_f_d2[5];
    double F_ave_f_d2[32];
    double F_lap_f_d2;
    double F_div_f_d2;
    double W_aveL_d3[5];
    double W_aveH_d3[5];
    double W_ave_f_d3[2048];
    double F_bar_f_d3[2048];    // Nearest power of 2 from 1715
    double W_f_d3[5];
    double F_ave_f_d3[32];
    double F_lap_f_d3;
    double F_div_f_d3;

// consToPrim1+deconvolve+consToPrim2+waveSpeedBound1+absMax+laplacian+increment+interpL_d1+interpH_d1+
// upwindState1+getFlux1+deconvolve_f_d1+getFlux2+smul_d1+lap_f_d1+inc_f_d1+div_f_d1+inc_rhs_d1+interpL_d2+
// interpH_d2+upwindState2+getFlux3+deconvolve_f_d2+getFlux4+smul_d2+lap_f_d2+inc_f_d2+div_f_d2+inc_rhs_d2+
// interpL_d3+interpH_d3+upwindState3+getFlux5+deconvolve_f_d3+getFlux6+lap_f_d3+inc_f_d3+div_f_d3+inc_rhs_d3+muldx

//    consToPrim1
#undef s0
#define s0(z,y,x) {\
W_bar(0,(z),(y),(x)) = U(0,(z),(y),(x));\
W_bar(1,(z),(y),(x)) = U(1,(z),(y),(x))/U(0,(z),(y),(x));\
W_bar(2,(z),(y),(x)) = U(2,(z),(y),(x))/U(0,(z),(y),(x));\
W_bar(3,(z),(y),(x)) = U(3,(z),(y),(x))/U(0,(z),(y),(x));\
W_bar(4,(z),(y),(x)) = (U(4,(z),(y),(x))-0.500000*U(0,(z),(y),(x))*(((U(1,(z),(y),(x))/U(0,(z),(y),(x)))*(U(1,(z),(y),(x))/U(0,(z),(y),(x))))+((U(2,(z),(y),(x))/U(0,(z),(y),(x)))*(U(2,(z),(y),(x))/U(0,(z),(y),(x))))+((U(3,(z),(y),(x))/U(0,(z),(y),(x)))*(U(3,(z),(y),(x))/U(0,(z),(y),(x))))))*(1.400000-1.000000);\
}
//    deconvolve
#undef s1
#define s1(c,z,y,x) u((c),(z),(y),(x))=(1.250000*U((c),(z),(y),(x)))+((-0.041667)*U((c),(z),(y),(x)+1))+((-0.041667)*U((c),(z),(y),(x)-1))+((-0.041667)*U((c),(z),(y)+1,(x)))+((-0.041667)*U((c),(z),(y)-1,(x)))+((-0.041667)*U((c),(z)+1,(y),(x)))+((-0.041667)*U((c),(z)-1,(y),(x)))
//    consToPrim2
#undef s2
#define s2(z,y,x) {\
W(0,(z),(y),(x)) = u(0,(z),(y),(x));\
W(1,(z),(y),(x)) = u(1,(z),(y),(x))/u(0,(z),(y),(x));\
W(2,(z),(y),(x)) = u(2,(z),(y),(x))/u(0,(z),(y),(x));\
W(3,(z),(y),(x)) = u(3,(z),(y),(x))/u(0,(z),(y),(x));\
W(4,(z),(y),(x)) = (u(4,(z),(y),(x))-0.500000*u(0,(z),(y),(x))*(((u(1,(z),(y),(x))/u(0,(z),(y),(x)))*(u(1,(z),(y),(x))/u(0,(z),(y),(x))))+((u(2,(z),(y),(x))/u(0,(z),(y),(x)))*(u(2,(z),(y),(x))/u(0,(z),(y),(x))))+((u(3,(z),(y),(x))/u(0,(z),(y),(x)))*(u(3,(z),(y),(x))/u(0,(z),(y),(x))))))*(1.400000-1.000000);\
}
// waveSpeedBound1
#undef s3
#define s3(z,y,x) {\
umax((z),(y),(x))=(3.000000*sqrt(1.400000*W(4,(z),(y),(x))/W(0,(z),(y),(x))))+W(1,(z),(y),(x))+W(2,(z),(y),(x))+W(3,(z),(y),(x));\
}
//     absMax
#undef s4
#define s4(z,y,x) absmax(retval,umax((z),(y),(x)))
//   laplacian
#undef s5
#define s5(c,z,y,x) W_ave((c),(z),(y),(x))=((-0.250000)*W_bar((c),(z),(y),(x)))+(0.041667*W_bar((c),(z),(y),(x)+1))+(0.041667*W_bar((c),(z),(y),(x)-1))+(0.041667*W_bar((c),(z),(y)+1,(x)))+(0.041667*W_bar((c),(z),(y)-1,(x)))+(0.041667*W_bar((c),(z)+1,(y),(x)))+(0.041667*W_bar((c),(z)-1,(y),(x)))
//    increment
#undef s6
#define s6(c,z,y,x) W_ave((c),(z),(y),(x))+=W((c),(z),(y),(x))
// interpL_d1
#undef s7
#define s7(c,z,y,x) W_aveL_d1((c),(z),(y),(x))=(0.033333*W_ave((c),(z),(y),(x)-3))+((-0.050000)*W_ave((c),(z),(y),(x)+1))+((-0.216667)*W_ave((c),(z),(y),(x)-2))+(0.450000*W_ave((c),(z),(y),(x)))+(0.783333*W_ave((c),(z),(y),(x)-1))
// interpH_d1
#undef s8
#define s8(c,z,y,x) W_aveH_d1((c),(z),(y),(x))=((-0.050000)*W_ave((c),(z),(y),(x)-2))+(0.033333*W_ave((c),(z),(y),(x)+2))+(0.450000*W_ave((c),(z),(y),(x)-1))+((-0.216667)*W_ave((c),(z),(y),(x)+1))+(0.783333*W_ave((c),(z),(y),(x)))
//    upwindState1
#undef s9
#define s9(z,y,x) {\
W_ave_f_d1(0,(z),(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d1(0,(z),(y),(x))) : (W_aveH_d1(0,(z),(y),(x)));\
W_ave_f_d1(1,(z),(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d1(1,(z),(y),(x))) : (W_aveH_d1(1,(z),(y),(x)));\
W_ave_f_d1(2,(z),(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d1(2,(z),(y),(x))) : (W_aveH_d1(2,(z),(y),(x)));\
W_ave_f_d1(3,(z),(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d1(3,(z),(y),(x))) : (W_aveH_d1(3,(z),(y),(x)));\
W_ave_f_d1(4,(z),(y),(x))=(0.000000 < ((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d1(4,(z),(y),(x))) : (W_aveH_d1(4,(z),(y),(x)));\
W_ave_f_d1(0,(z),(y),(x))+=(0.000000 < sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000) ? (((((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000+((((((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))*(W_aveL_d1(0+1,(z),(y),(x))-W_aveH_d1(0+1,(z),(y),(x))))*0.500000-W_ave_f_d1(4,(z),(y),(x))))/((sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000)))*sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000)))) : (0.000000);\
W_ave_f_d1(0+1,(z),(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000) ? (((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))) : (W_ave_f_d1(0+1,(z),(y),(x)));\
W_ave_f_d1(4,(z),(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000+((W_aveL_d1(4,(z),(y),(x))-W_aveH_d1(4,(z),(y),(x))))/((2.000000*(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d1(0+1,(z),(y),(x))+W_aveH_d1(0+1,(z),(y),(x))))*0.500000) ? (((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000+((((((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d1(4,(z),(y),(x))+W_aveH_d1(4,(z),(y),(x))))*0.500000))/(((W_aveL_d1(0,(z),(y),(x))+W_aveH_d1(0,(z),(y),(x))))*0.500000))))*(W_aveL_d1(0+1,(z),(y),(x))-W_aveH_d1(0+1,(z),(y),(x))))*0.500000) : (W_ave_f_d1(4,(z),(y),(x)));\
}
//    getFlux1
#undef s10
#define s10(z,y,x) {\
F_bar_f_d1(0,(z),(y),(x))=(W_ave_f_d1(0+1,(z),(y),(x)))*W_ave_f_d1(0,(z),(y),(x));\
F_bar_f_d1(1,(z),(y),(x))=W_ave_f_d1(1,(z),(y),(x))*F_bar_f_d1(0,(z),(y),(x));\
F_bar_f_d1(2,(z),(y),(x))=W_ave_f_d1(2,(z),(y),(x))*F_bar_f_d1(0,(z),(y),(x));\
F_bar_f_d1(3,(z),(y),(x))=W_ave_f_d1(3,(z),(y),(x))*F_bar_f_d1(0,(z),(y),(x));\
F_bar_f_d1(0+1,(z),(y),(x))+=W_ave_f_d1(4,(z),(y),(x));\
F_bar_f_d1(4,(z),(y),(x))=((1.400000/(1.400000-1))*W_ave_f_d1(0+1,(z),(y),(x)))*W_ave_f_d1(4,(z),(y),(x))+0.500000*F_bar_f_d1(0,(z),(y),(x))*((W_ave_f_d1(1,(z),(y),(x))*W_ave_f_d1(1,(z),(y),(x)))+(W_ave_f_d1(2,(z),(y),(x))*W_ave_f_d1(2,(z),(y),(x)))+(W_ave_f_d1(3,(z),(y),(x))*W_ave_f_d1(3,(z),(y),(x))));\
}
// deconvolve_f_d1
#undef s11
#define s11(c,z,y,x) W_f_d1((c),(z),(y),(x))=(1.166667*W_ave_f_d1((c),(z),(y),(x)))+((-0.041667)*W_ave_f_d1((c),(z),(y)+1,(x)))+((-0.041667)*W_ave_f_d1((c),(z),(y)-1,(x)))+((-0.041667)*W_ave_f_d1((c),(z)+1,(y),(x)))+((-0.041667)*W_ave_f_d1((c),(z)-1,(y),(x)))
// getFlux2
#undef s12
#define s12(z,y,x) {\
F_ave_f_d1(0,(z),(y),(x))=(W_f_d1(0+1,(z),(y),(x)))*W_f_d1(0,(z),(y),(x));\
F_ave_f_d1(1,(z),(y),(x))=W_f_d1(1,(z),(y),(x))*F_ave_f_d1(0,(z),(y),(x));\
F_ave_f_d1(2,(z),(y),(x))=W_f_d1(2,(z),(y),(x))*F_ave_f_d1(0,(z),(y),(x));\
F_ave_f_d1(3,(z),(y),(x))=W_f_d1(3,(z),(y),(x))*F_ave_f_d1(0,(z),(y),(x));\
F_ave_f_d1(0+1,(z),(y),(x))+=W_f_d1(4,(z),(y),(x));\
F_ave_f_d1(4,(z),(y),(x))=((1.400000/(1.400000-1))*W_f_d1(0+1,(z),(y),(x)))*W_f_d1(4,(z),(y),(x))+0.500000*F_ave_f_d1(0,(z),(y),(x))*((W_f_d1(1,(z),(y),(x))*W_f_d1(1,(z),(y),(x)))+(W_f_d1(2,(z),(y),(x))*W_f_d1(2,(z),(y),(x)))+(W_f_d1(3,(z),(y),(x))*W_f_d1(3,(z),(y),(x))));\
}
//    smul_d1
#undef s13
#define s13(c,z,y,x) F_bar_f_d1((c),(z),(y),(x))*=0.041667
//   lap_f_d1
#undef s14
#define s14(c,z,y,x) F_lap_f_d1((c),(z),(y),(x))=((-0.166667)*F_bar_f_d1((c),(z),(y),(x)))+(0.041667*F_bar_f_d1((c),(z),(y)+1,(x)))+(0.041667*F_bar_f_d1((c),(z),(y)-1,(x)))+(0.041667*F_bar_f_d1((c),(z)+1,(y),(x)))+(0.041667*F_bar_f_d1((c),(z)-1,(y),(x)))
// inc_f_d1
#undef s15
#define s15(c,z,y,x) F_ave_f_d1((c),(z),(y),(x))+=F_lap_f_d1((c),(z),(y),(x))
//    div_f_d1
#undef s16
#define s16(c,z,y,x) F_div_f_d1((c),(z),(y),(x))=((-1.000000)*F_ave_f_d1((c),(z),(y),(x)))+(1.000000*F_ave_f_d1((c),(z),(y),(x)+1))
//    inc_rhs_d1
#undef s17
#define s17(c,z,y,x) rhs((c),(z),(y),(x))+=F_div_f_d1((c),(z),(y),(x))
// interpL_d2
#undef s18
#define s18(c,z,y,x) W_aveL_d2((c),(z),(y),(x))=(0.033333*W_ave((c),(z),(y)-3,(x)))+((-0.050000)*W_ave((c),(z),(y)+1,(x)))+((-0.216667)*W_ave((c),(z),(y)-2,(x)))+(0.450000*W_ave((c),(z),(y),(x)))+(0.783333*W_ave((c),(z),(y)-1,(x)))
//   interpH_d2
#undef s19
#define s19(c,z,y,x) W_aveH_d2((c),(z),(y),(x))=((-0.050000)*W_ave((c),(z),(y)-2,(x)))+(0.033333*W_ave((c),(z),(y)+2,(x)))+(0.450000*W_ave((c),(z),(y)-1,(x)))+((-0.216667)*W_ave((c),(z),(y)+1,(x)))+(0.783333*W_ave((c),(z),(y),(x)))
// upwindState2
#undef s20
#define s20(z,y,x) {\
W_ave_f_d2(0,(z),(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d2(0,(z),(y),(x))) : (W_aveH_d2(0,(z),(y),(x)));\
W_ave_f_d2(1,(z),(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d2(1,(z),(y),(x))) : (W_aveH_d2(1,(z),(y),(x)));\
W_ave_f_d2(2,(z),(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d2(2,(z),(y),(x))) : (W_aveH_d2(2,(z),(y),(x)));\
W_ave_f_d2(3,(z),(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d2(3,(z),(y),(x))) : (W_aveH_d2(3,(z),(y),(x)));\
W_ave_f_d2(4,(z),(y),(x))=(0.000000 < ((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d2(4,(z),(y),(x))) : (W_aveH_d2(4,(z),(y),(x)));\
W_ave_f_d2(0,(z),(y),(x))+=(0.000000 < sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000) ? (((((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000+((((((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))*(W_aveL_d2(1+1,(z),(y),(x))-W_aveH_d2(1+1,(z),(y),(x))))*0.500000-W_ave_f_d2(4,(z),(y),(x))))/((sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000)))*sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000)))) : (0.000000);\
W_ave_f_d2(1+1,(z),(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000) ? (((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))) : (W_ave_f_d2(1+1,(z),(y),(x)));\
W_ave_f_d2(4,(z),(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000+((W_aveL_d2(4,(z),(y),(x))-W_aveH_d2(4,(z),(y),(x))))/((2.000000*(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d2(1+1,(z),(y),(x))+W_aveH_d2(1+1,(z),(y),(x))))*0.500000) ? (((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000+((((((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d2(4,(z),(y),(x))+W_aveH_d2(4,(z),(y),(x))))*0.500000))/(((W_aveL_d2(0,(z),(y),(x))+W_aveH_d2(0,(z),(y),(x))))*0.500000))))*(W_aveL_d2(1+1,(z),(y),(x))-W_aveH_d2(1+1,(z),(y),(x))))*0.500000) : (W_ave_f_d2(4,(z),(y),(x)));\
}
// getFlux3
#undef s21
#define s21(z,y,x) {\
F_bar_f_d2(0,(z),(y),(x))=(W_ave_f_d2(1+1,(z),(y),(x)))*W_ave_f_d2(0,(z),(y),(x));\
F_bar_f_d2(1,(z),(y),(x))=W_ave_f_d2(1,(z),(y),(x))*F_bar_f_d2(0,(z),(y),(x));\
F_bar_f_d2(2,(z),(y),(x))=W_ave_f_d2(2,(z),(y),(x))*F_bar_f_d2(0,(z),(y),(x));\
F_bar_f_d2(3,(z),(y),(x))=W_ave_f_d2(3,(z),(y),(x))*F_bar_f_d2(0,(z),(y),(x));\
F_bar_f_d2(1+1,(z),(y),(x))+=W_ave_f_d2(4,(z),(y),(x));\
F_bar_f_d2(4,(z),(y),(x))=((1.400000/(1.400000-1))*W_ave_f_d2(1+1,(z),(y),(x)))*W_ave_f_d2(4,(z),(y),(x))+0.500000*F_bar_f_d2(0,(z),(y),(x))*((W_ave_f_d2(1,(z),(y),(x))*W_ave_f_d2(1,(z),(y),(x)))+(W_ave_f_d2(2,(z),(y),(x))*W_ave_f_d2(2,(z),(y),(x)))+(W_ave_f_d2(3,(z),(y),(x))*W_ave_f_d2(3,(z),(y),(x))));\
}
//    deconvolve_f_d2
#undef s22
#define s22(c,z,y,x) W_f_d2((c),(z),(y),(x))=(1.166667*W_ave_f_d2((c),(z),(y),(x)))+((-0.041667)*W_ave_f_d2((c),(z),(y),(x)+1))+((-0.041667)*W_ave_f_d2((c),(z),(y),(x)-1))+((-0.041667)*W_ave_f_d2((c),(z)+1,(y),(x)))+((-0.041667)*W_ave_f_d2((c),(z)-1,(y),(x)))
// getFlux4
#undef s23
#define s23(z,y,x) {\
F_ave_f_d2(0,(z),(y),(x))=(W_f_d2(1+1,(z),(y),(x)))*W_f_d2(0,(z),(y),(x));\
F_ave_f_d2(1,(z),(y),(x))=W_f_d2(1,(z),(y),(x))*F_ave_f_d2(0,(z),(y),(x));\
F_ave_f_d2(2,(z),(y),(x))=W_f_d2(2,(z),(y),(x))*F_ave_f_d2(0,(z),(y),(x));\
F_ave_f_d2(3,(z),(y),(x))=W_f_d2(3,(z),(y),(x))*F_ave_f_d2(0,(z),(y),(x));\
F_ave_f_d2(1+1,(z),(y),(x))+=W_f_d2(4,(z),(y),(x));\
F_ave_f_d2(4,(z),(y),(x))=((1.400000/(1.400000-1))*W_f_d2(1+1,(z),(y),(x)))*W_f_d2(4,(z),(y),(x))+0.500000*F_ave_f_d2(0,(z),(y),(x))*((W_f_d2(1,(z),(y),(x))*W_f_d2(1,(z),(y),(x)))+(W_f_d2(2,(z),(y),(x))*W_f_d2(2,(z),(y),(x)))+(W_f_d2(3,(z),(y),(x))*W_f_d2(3,(z),(y),(x))));\
}
// smul_d2
#undef s24
#define s24(c,z,y,x) F_bar_f_d2((c),(z),(y),(x))*=0.041667
// lap_f_d2
#undef s25
#define s25(c,z,y,x) F_lap_f_d2((c),(z),(y),(x))=((-0.166667)*F_bar_f_d2((c),(z),(y),(x)))+(0.041667*F_bar_f_d2((c),(z),(y),(x)+1))+(0.041667*F_bar_f_d2((c),(z),(y),(x)-1))+(0.041667*F_bar_f_d2((c),(z)+1,(y),(x)))+(0.041667*F_bar_f_d2((c),(z)-1,(y),(x)))
// inc_f_d2
#undef s26
#define s26(c,z,y,x) F_ave_f_d2((c),(z),(y),(x))+=F_lap_f_d2((c),(z),(y),(x))
// div_f_d2
#undef s27
#define s27(c,z,y,x) F_div_f_d2((c),(z),(y),(x))=((-1.000000)*F_ave_f_d2((c),(z),(y),(x)))+(1.000000*F_ave_f_d2((c),(z),(y)+1,(x)))
// inc_rhs_d2
#undef s28
#define s28(c,z,y,x) rhs((c),(z),(y),(x))+=F_div_f_d2((c),(z),(y),(x))
// interpL_d3
#undef s29
#define s29(c,z,y,x) W_aveL_d3((c),(z),(y),(x))=(0.033333*W_ave((c),(z)-3,(y),(x)))+((-0.050000)*W_ave((c),(z)+1,(y),(x)))+((-0.216667)*W_ave((c),(z)-2,(y),(x)))+(0.450000*W_ave((c),(z),(y),(x)))+(0.783333*W_ave((c),(z)-1,(y),(x)))
// interpH_d3
#undef s30
#define s30(c,z,y,x) W_aveH_d3((c),(z),(y),(x))=((-0.050000)*W_ave((c),(z)-2,(y),(x)))+(0.033333*W_ave((c),(z)+2,(y),(x)))+(0.450000*W_ave((c),(z)-1,(y),(x)))+((-0.216667)*W_ave((c),(z)+1,(y),(x)))+(0.783333*W_ave((c),(z),(y),(x)))
// upwindState3
#undef s31
#define s31(z,y,x) {\
W_ave_f_d3(0,(z),(y),(x))=(0.000000 < ((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d3(0,(z),(y),(x))) : (W_aveH_d3(0,(z),(y),(x)));\
W_ave_f_d3(1,(z),(y),(x))=(0.000000 < ((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d3(1,(z),(y),(x))) : (W_aveH_d3(1,(z),(y),(x)));\
W_ave_f_d3(2,(z),(y),(x))=(0.000000 < ((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d3(2,(z),(y),(x))) : (W_aveH_d3(2,(z),(y),(x)));\
W_ave_f_d3(3,(z),(y),(x))=(0.000000 < ((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d3(3,(z),(y),(x))) : (W_aveH_d3(3,(z),(y),(x)));\
W_ave_f_d3(4,(z),(y),(x))=(0.000000 < ((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))) ? (W_aveL_d3(4,(z),(y),(x))) : (W_aveH_d3(4,(z),(y),(x)));\
W_ave_f_d3(0,(z),(y),(x))+=(0.000000 < sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000) ? (((((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000+((((((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))*(W_aveL_d3(2+1,(z),(y),(x))-W_aveH_d3(2+1,(z),(y),(x))))*0.500000-W_ave_f_d3(4,(z),(y),(x))))/((sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000)))*sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000)))) : (0.000000);\
W_ave_f_d3(2+1,(z),(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000) ? (((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))) : (W_ave_f_d3(2+1,(z),(y),(x)));\
W_ave_f_d3(4,(z),(y),(x))=(0.000000 < sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))+(-sgn(((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000+((W_aveL_d3(4,(z),(y),(x))-W_aveH_d3(4,(z),(y),(x))))/((2.000000*(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))))*((W_aveL_d3(2+1,(z),(y),(x))+W_aveH_d3(2+1,(z),(y),(x))))*0.500000) ? (((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000+((((((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))*(sqrt((1.400000*(((W_aveL_d3(4,(z),(y),(x))+W_aveH_d3(4,(z),(y),(x))))*0.500000))/(((W_aveL_d3(0,(z),(y),(x))+W_aveH_d3(0,(z),(y),(x))))*0.500000))))*(W_aveL_d3(2+1,(z),(y),(x))-W_aveH_d3(2+1,(z),(y),(x))))*0.500000) : (W_ave_f_d3(4,(z),(y),(x)));\
}
// getFlux5
#undef s32
#define s32(z,y,x) {\
F_bar_f_d3(0,(z),(y),(x))=(W_ave_f_d3(2+1,(z),(y),(x)))*W_ave_f_d3(0,(z),(y),(x));\
F_bar_f_d3(1,(z),(y),(x))=W_ave_f_d3(1,(z),(y),(x))*F_bar_f_d3(0,(z),(y),(x));\
F_bar_f_d3(2,(z),(y),(x))=W_ave_f_d3(2,(z),(y),(x))*F_bar_f_d3(0,(z),(y),(x));\
F_bar_f_d3(3,(z),(y),(x))=W_ave_f_d3(3,(z),(y),(x))*F_bar_f_d3(0,(z),(y),(x));\
F_bar_f_d3(2+1,(z),(y),(x))+=W_ave_f_d3(4,(z),(y),(x));\
F_bar_f_d3(4,(z),(y),(x))=((1.400000/(1.400000-1))*W_ave_f_d3(2+1,(z),(y),(x)))*W_ave_f_d3(4,(z),(y),(x))+0.500000*F_bar_f_d3(0,(z),(y),(x))*((W_ave_f_d3(1,(z),(y),(x))*W_ave_f_d3(1,(z),(y),(x)))+(W_ave_f_d3(2,(z),(y),(x))*W_ave_f_d3(2,(z),(y),(x)))+(W_ave_f_d3(3,(z),(y),(x))*W_ave_f_d3(3,(z),(y),(x))));\
}
// deconvolve_f_d3
#undef s33
#define s33(c,z,y,x) W_f_d3((c),(z),(y),(x))=(1.166667*W_ave_f_d3((c),(z),(y),(x)))+((-0.041667)*W_ave_f_d3((c),(z),(y),(x)+1))+((-0.041667)*W_ave_f_d3((c),(z),(y),(x)-1))+((-0.041667)*W_ave_f_d3((c),(z),(y)+1,(x)))+((-0.041667)*W_ave_f_d3((c),(z),(y)-1,(x)))
// getFlux6
#undef s34
#define s34(z,y,x) {\
F_ave_f_d3(0,(z),(y),(x))=(W_f_d3(2+1,(z),(y),(x)))*W_f_d3(0,(z),(y),(x));\
F_ave_f_d3(1,(z),(y),(x))=W_f_d3(1,(z),(y),(x))*F_ave_f_d3(0,(z),(y),(x));\
F_ave_f_d3(2,(z),(y),(x))=W_f_d3(2,(z),(y),(x))*F_ave_f_d3(0,(z),(y),(x));\
F_ave_f_d3(3,(z),(y),(x))=W_f_d3(3,(z),(y),(x))*F_ave_f_d3(0,(z),(y),(x));\
F_ave_f_d3(2+1,(z),(y),(x))+=W_f_d3(4,(z),(y),(x));\
F_ave_f_d3(4,(z),(y),(x))=((1.400000/(1.400000-1))*W_f_d3(2+1,(z),(y),(x)))*W_f_d3(4,(z),(y),(x))+0.500000*F_ave_f_d3(0,(z),(y),(x))*((W_f_d3(1,(z),(y),(x))*W_f_d3(1,(z),(y),(x)))+(W_f_d3(2,(z),(y),(x))*W_f_d3(2,(z),(y),(x)))+(W_f_d3(3,(z),(y),(x))*W_f_d3(3,(z),(y),(x))));\
}
// smul_d3
#undef s35
#define s35(c,z,y,x) F_bar_f_d3((c),(z),(y),(x))*=0.041667
// lap_f_d3
#undef s36
#define s36(c,z,y,x) F_lap_f_d3((c),(z),(y),(x))=((-0.166667)*F_bar_f_d3((c),(z),(y),(x)))+(0.041667*F_bar_f_d3((c),(z),(y),(x)+1))+(0.041667*F_bar_f_d3((c),(z),(y),(x)-1))+(0.041667*F_bar_f_d3((c),(z),(y)+1,(x)))+(0.041667*F_bar_f_d3((c),(z),(y)-1,(x)))
// inc_f_d3
#undef s37
#define s37(c,z,y,x) F_ave_f_d3((c),(z),(y),(x))+=F_lap_f_d3((c),(z),(y),(x))
// div_f_d3
#undef s38
#define s38(c,z,y,x) F_div_f_d3((c),(z),(y),(x))=((-1.000000)*F_ave_f_d3((c),(z),(y),(x)))+(1.000000*F_ave_f_d3((c),(z)+1,(y),(x)))
// inc_rhs_d3
#undef s39
#define s39(c,z,y,x) rhs((c),(z),(y),(x))+=F_div_f_d3((c),(z),(y),(x))
// muldx
#undef s40
#define s40(c,z,y,x) rhs((c),(z),(y),(x))*=-1.000000

for(t1 = -4; t1 <= 67; t1++) {
  for(t2 = -4; t2 <= 67; t2++) {
#pragma omp simd
    for(t3 = -4; t3 <= 67; t3++) {
      //z = t1, y = t2, x = t3;
        // Precompute future spatial dimensions.
        s0(t1,t2,t3);
        s0(t1,t2,t3+1);
        s0(t1,t2,t3+2);
        s0(t1,t2,t3+3);
        s0(t1,t2+1,t3);
        s0(t1,t2+2,t3);
        s0(t1,t2+3,t3);
        s0(t1+1,t2,t3);
        s0(t1+2,t2,t3);
        s0(t1+3,t2,t3);

      s1(0,t1,t2,t3);
      s1(1,t1,t2,t3);
      s1(2,t1,t2,t3);
      s1(3,t1,t2,t3);
      s1(4,t1,t2,t3);

      s2(t1,t2,t3);
      s3(t1,t2,t3);
      s4(t1,t2,t3);

        s5(0,t1,t2,t3);
        s5(1,t1,t2,t3);
        s5(2,t1,t2,t3);
        s5(3,t1,t2,t3);
        s5(4,t1,t2,t3);

        s6(0,t1,t2,t3);
        s6(1,t1,t2,t3);
        s6(2,t1,t2,t3);
        s6(3,t1,t2,t3);
        s6(4,t1,t2,t3);

        s7(0,t1,t2,t3);
        s7(1,t1,t2,t3);
        s7(2,t1,t2,t3);
        s7(3,t1,t2,t3);
        s7(4,t1,t2,t3);

        s8(0,t1,t2,t3);
        s8(1,t1,t2,t3);
        s8(2,t1,t2,t3);
        s8(3,t1,t2,t3);
        s8(4,t1,t2,t3);

      s9(t1,t2,t3);
      s10(t1,t2,t3);

        s11(0,t1,t2,t3);
        s11(1,t1,t2,t3);
        s11(2,t1,t2,t3);
        s11(3,t1,t2,t3);
        s11(4,t1,t2,t3);

        s13(0,t1,t2,t3);
        s13(1,t1,t2,t3);
        s13(2,t1,t2,t3);
        s13(3,t1,t2,t3);
        s13(4,t1,t2,t3);

      s12(t1,t2,t3);

        s14(0,t1,t2,t3);
        s14(1,t1,t2,t3);
        s14(2,t1,t2,t3);
        s14(3,t1,t2,t3);
        s14(4,t1,t2,t3);

        s15(0,t1,t2,t3);
        s15(1,t1,t2,t3);
        s15(2,t1,t2,t3);
        s15(3,t1,t2,t3);
        s15(4,t1,t2,t3);

        s16(0,t1,t2,t3);
        s16(1,t1,t2,t3);
        s16(2,t1,t2,t3);
        s16(3,t1,t2,t3);
        s16(4,t1,t2,t3);

        // TODO: Output (rhs) cannot be resized so this guard is necessary, unless I use a temp buffer and copy at end.
        if ((t1 >= 0 && t1 <= 63) && (t2 >= 0 && t2 <= 63) && (t3 >= 0 && t3 <= 63)) {
            s17(0,t1,t2,t3);
            s17(1,t1,t2,t3);
            s17(2,t1,t2,t3);
            s17(3,t1,t2,t3);
            s17(4,t1,t2,t3);
        }

        s18(0,t1,t2,t3);
        s18(1,t1,t2,t3);
        s18(2,t1,t2,t3);
        s18(3,t1,t2,t3);
        s18(4,t1,t2,t3);

        s19(0,t1,t2,t3);
        s19(1,t1,t2,t3);
        s19(2,t1,t2,t3);
        s19(3,t1,t2,t3);
        s19(4,t1,t2,t3);

      s20(t1,t2,t3);
      s21(t1,t2,t3);

        s22(0,t1,t2,t3);
        s22(1,t1,t2,t3);
        s22(2,t1,t2,t3);
        s22(3,t1,t2,t3);
        s22(4,t1,t2,t3);

        s23(t1,t2,t3);

        s24(0,t1,t2,t3);
        s24(1,t1,t2,t3);
        s24(2,t1,t2,t3);
        s24(3,t1,t2,t3);
        s24(4,t1,t2,t3);

        s25(0,t1,t2,t3);
        s25(1,t1,t2,t3);
        s25(2,t1,t2,t3);
        s25(3,t1,t2,t3);
        s25(4,t1,t2,t3);

        s26(0,t1,t2,t3);
        s26(1,t1,t2,t3);
        s26(2,t1,t2,t3);
        s26(3,t1,t2,t3);
        s26(4,t1,t2,t3);

        s27(0,t1,t2,t3);
        s27(1,t1,t2,t3);
        s27(2,t1,t2,t3);
        s27(3,t1,t2,t3);
        s27(4,t1,t2,t3);

        // TODO: Output (rhs) cannot be resized so this guard is necessary, unless I use a temp buffer and copy at end.
        if ((t1 >= 0 && t1 <= 63) && (t2 >= 0 && t2 <= 63) && (t3 >= 0 && t3 <= 63)) {
            s28(0,t1,t2,t3);
            s28(1,t1,t2,t3);
            s28(2,t1,t2,t3);
            s28(3,t1,t2,t3);
            s28(4,t1,t2,t3);
        }

        s29(0,t1,t2,t3);
        s29(1,t1,t2,t3);
        s29(2,t1,t2,t3);
        s29(3,t1,t2,t3);
        s29(4,t1,t2,t3);

        s30(0,t1,t2,t3);
        s30(1,t1,t2,t3);
        s30(2,t1,t2,t3);
        s30(3,t1,t2,t3);
        s30(4,t1,t2,t3);

      s31(t1,t2,t3);
      s32(t1,t2,t3);

        s33(0,t1,t2,t3);
        s33(1,t1,t2,t3);
        s33(2,t1,t2,t3);
        s33(3,t1,t2,t3);
        s33(4,t1,t2,t3);

        s34(t1,t2,t3);

        s35(0,t1,t2,t3);
        s35(1,t1,t2,t3);
        s35(2,t1,t2,t3);
        s35(3,t1,t2,t3);
        s35(4,t1,t2,t3);

        s36(0,t1,t2,t3);
        s36(1,t1,t2,t3);
        s36(2,t1,t2,t3);
        s36(3,t1,t2,t3);
        s36(4,t1,t2,t3);

        s37(0,t1,t2,t3);
        s37(1,t1,t2,t3);
        s37(2,t1,t2,t3);
        s37(3,t1,t2,t3);
        s37(4,t1,t2,t3);

        s38(0,t1,t2,t3);
        s38(1,t1,t2,t3);
        s38(2,t1,t2,t3);
        s38(3,t1,t2,t3);
        s38(4,t1,t2,t3);

        // TODO: Output (rhs) cannot be resized so this guard is necessary, unless I use a temp buffer and copy at end.
        if ((t1 >= 0 && t1 <= 63) && (t2 >= 0 && t2 <= 63) && (t3 >= 0 && t3 <= 63)) {
            s39(0,t1,t2,t3);
            s39(1,t1,t2,t3);
            s39(2,t1,t2,t3);
            s39(3,t1,t2,t3);
            s39(4,t1,t2,t3);

            s40(0,t1,t2,t3);
            s40(1,t1,t2,t3);
            s40(2,t1,t2,t3);
            s40(3,t1,t2,t3);
            s40(4,t1,t2,t3);
        }
    }
  }
}

//    free(U);
//    fprintf(stderr, "retval: %g\n", retval);

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
#undef U
#undef W_bar
#undef u
#undef W
#undef umax
#undef W_ave
#undef W_aveL_d1
#undef W_aveH_d1
#undef W_ave_f_d1
#undef F_bar_f_d1
#undef W_f_d1
#undef F_ave_f_d1
#undef F_lap_f_d1
#undef F_div_f_d1
#undef rhs
#undef W_aveL_d2
#undef W_aveH_d2
#undef W_ave_f_d2
#undef F_bar_f_d2
#undef W_f_d2
#undef F_ave_f_d2
#undef F_lap_f_d2
#undef F_div_f_d2
#undef W_aveL_d3
#undef W_aveH_d3
#undef W_ave_f_d3
#undef F_bar_f_d3
#undef W_f_d3
#undef F_ave_f_d3
#undef F_lap_f_d3
#undef F_div_f_d3
