#ifndef _TESTMAPPING_H_
#define _TESTMAPPING_H_
#include "Proto.H"
using namespace std;
using namespace Proto;
class PatchMap
{
public:  
  double m_h;
  array<array<double,DIM>,DIM> m_array;
  PatchMap(const double& a_h)
  {
    m_h = a_h;
    for (int dir0 =0;dir0 < DIM;dir0++)
      for (int dir1=0;dir1 < DIM;dir1++)
        {
          m_array[dir0][dir1]=Point::Basis(dir0,1)[dir1];
        }
  };
  inline
  BoxData<double,DIM> map(const Box& a_bx,int a_nGhost)
  {
    Box bxmap =a_bx.grow(a_nGhost);
    for (int dir = 0;dir <DIM;dir++)
      {
        bxmap =bxmap.extrude(dir);
      }
    BoxData<double,DIM> ret(bxmap);
    forallInPlace_p(
                    [] PROTO_LAMBDA(const Point& a_pt,
                                    Var<double,DIM>& a_pos,
                                    const array<array<double,DIM>,DIM>& a_array,
                                    const double& a_h)
                    {
                      for (int dir0 = 0; dir0 < DIM; dir0++)
                        for (int dir1 = 0; dir1 < DIM; dir1++)
                        {
                          a_pos(dir0) += a_pt[dir1]*a_array[dir0][dir1]*a_h;
                        }
                    },ret,m_array,m_h);
    return ret;
  };
};
BoxData<double,DIM> fAv(
                        const BoxData<double,DIM>& a_X,
                        const Point& a_k,
                        int a_dir)
{
  HDF5Handler h5;
  auto fluxCorner = forall<double,DIM>
    (
     [ ] PROTO_LAMBDA (
                       Var<double,DIM>& a_flux,                       
                       const Var<double,DIM>& a_pos,
                       const Point& a_k)
     { 
       for (int dir1 = 0; dir1 < DIM; dir1++)
         {
           a_flux(dir1) = a_pos(dir1);
           /* a_flux(dir1) = 1.0;
                     for (int dir2 = 0; dir2 < DIM; dir2++)
             {
               if (dir2==dir1)
                 {
                   a_flux(dir1) *= 2*M_PI*a_k[dir1]
                     *cos(2*M_PI*a_k[dir1]*a_pos(dir1));
                 }
               else
                 {
                 a_flux(dir1) *= sin(2*M_PI*a_k[dir2]*a_pos(dir2)); 
                   //cout << "dir1,dir2,flux(dir1) "<< dir1<< " , "<<dir2 << " , " << sin(2*M_PI*a_k[dir2]*a_pos(dir2)) << endl;
                 }
                 //cout << "dir1,dir2,flux(dir1) "<< dir1<< " , "<<dir2 << " , " << a_flux(dir1) << endl; 
                 } */
         }
     },a_X,a_k);
  string strdir = "fluxCornerDir" + to_string(a_dir);
  h5.writePatch(1.0,fluxCorner,strdir);
  BoxData<double,DIM> ret = Stencil<double>::cornersToFaces(a_dir,4)(fluxCorner);
  strdir = "retFluxDir" + to_string(a_dir);
  h5.writePatch(1.0,ret,strdir);
  return ret;
}
    BoxData<double> divFExact(const Box& a_bx,
                              const BoxData<double,DIM>& a_X,
                              const Point& a_k)
{
  BoxData<double> ret(a_bx);
  ret.setToZero();
  auto lapCorner = forall<double>
    (
     [ ] PROTO_LAMBDA (
                       Var<double>& a_ret,
                       const Var<double,DIM >& a_X,
                       const Point& a_k)
     {
       int ksq = 0;
       a_ret(0) = DIM*1.0;
       /*for (int dir = 0; dir < DIM;dir++)
         {
           ksq+= a_k[dir]*a_k[dir];
           a_ret(0) *= sin(2*M_PI*a_X(dir));
         }
       a_ret(0)*= -4.0*M_PI*M_PI*ksq;*/
     },a_X,a_k);
  ret += Stencil<double>::cornersToCells(4)(lapCorner);
  return ret;
}
#endif
