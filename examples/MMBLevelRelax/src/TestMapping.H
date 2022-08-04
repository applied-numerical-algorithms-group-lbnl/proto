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
BoxData<double,DIM> fAv(const Box& a_bx,
                        const BoxData<double,DIM>& a_X,
                        int a_dir)
{
  auto ret = forall<double,DIM >
    (
     [ ] PROTO_LAMBDA (
                       Var<double,DIM>& a_ret,
                       const Var<double,DIM>& a_pos,
                       int a_dir)
     {
       for (int dir = 0; dir < DIM; dir++)
         {
           if (dir == a_dir)
             {
               a_ret(dir) = (1.0/DIM)*a_pos(dir);
             }
           else
             {
               a_ret(dir) = 0.;
             }
         }
     },a_X,a_dir);
  return ret;
}
BoxData<double> divFExact(const Box& a_bx, const BoxData<double,DIM>& a_X)
{
  BoxData<double> ret(a_bx);
  forallInPlace(
                [ ] PROTO_LAMBDA (
                                  Var<double>& a_ret,
                                  const Var<double,DIM >& a_X)
                {
                  a_ret(0) = 1.0;
                },ret,a_X);
    return ret;                                   
}
#endif
