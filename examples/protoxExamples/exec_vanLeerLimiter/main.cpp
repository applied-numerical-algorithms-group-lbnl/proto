#include "Proto.H"
#define NUMCOMPS 1
#define NGHOST 2
#define MAXLEV 2
#define ORDER 4
using namespace Proto;
// define van Leer pointwise function.

template<typename T,unsigned int C, MemType MEM>
PROTO_KERNEL_START
void vanleer_F(
               Var<T,C,MEM>& a_dwlim,
               Var<T,C,MEM>& a_dwlo,
               Var<T,C,MEM>& a_dwhi,
               T a_h)
{
  for (int comp = 0; comp < C;comp++)
    {
      T dwlo = a_dwlo(comp);
      T dwhi = a_dwhi(comp);
      T dwlim = min(
        min(2*abs(dwlo),2*abs(dwhi)),abs(dwlo + dwhi)*.5);
      //int slo,shi;
      //if (std::signbit(dwlo)) slo = -1; else slo = 1;
      //if (std::signbit(dwhi)) shi = -1; else shi = 1;
      if (dwlo*dwhi <= 0.)
        {
          dwlim = 0.;
        }
      else if (dwlo < 0)
        {
          dwlim *= -1.0;
        }
      a_dwlim(comp) = dwlim/a_h;
    }
}
PROTO_KERNEL_END(vanleer_F,vanleer)
int main(int argc, char *argv[])
{

  int nx;
  int dir;
  cout << "input nx, dir" << endl;
  cin >> nx >> dir;
  PR_TIMER_SETFILE("vanLeerOperator" + to_string(nx)+"_"+to_string(dir)+".time.table");
  PR_TIMERS("Laplacian_stencil");
  Array<double,MAXLEV> errorMax,errorL1;
  for (int lev = 0; lev < MAXLEV;lev++)
    {
      // initialize data.
      Box bx(Point::Ones(-NGHOST),Point::Ones(nx+NGHOST));
      Box bx0(Point::Ones(1),Point::Ones(nx-1));
      BoxData<double,NUMCOMPS,HOST> u(bx),LOfuE(bx);
      double h = 1.0/nx;
      forallInPlace_p([ ]PROTO_LAMBDA(Point a_pt,
                                      Var<double,NUMCOMPS>& a_u,
                                      Var<double,NUMCOMPS>& a_lofue,
                                      double a_h,
                                      int a_dir)
                      {
                        for (int comp = 0; comp < NUMCOMPS;comp++)
                          {
                            a_u(comp) = 1.;
                            a_lofue(comp) = 1.;
                            for (int dir = 0; dir < DIM; dir++)
                              {
                                a_u(comp) *= sin(a_pt[dir]*M_PI*a_h + .25*M_PI*a_h);
                                if (dir == a_dir)
                                  {
                                    a_lofue(comp) *= M_PI*cos(a_pt[dir]*M_PI*a_h + .25*M_PI*a_h);
                                  }
                                else
                                  {
                                    a_lofue(comp) *= sin(a_pt[dir]*M_PI*a_h  + .25*M_PI*a_h);

                                  }
                              }
                          }
                      },bx,u,LOfuE,h,dir);
      
      // apply operator.
      BoxData<double,NUMCOMPS,HOST> LOfu;
      {
        PR_TIMERS("timed section");
        BoxData<double,NUMCOMPS,HOST> du =
          ((1.0)*Shift(Point::Basis(dir,1)) + (-1.0)*Shift(Point::Zeros()))(u);
        BoxData<double,NUMCOMPS> dulo = alias(du,Point::Basis(dir,1));
        LOfu = forall<double,NUMCOMPS,HOST>(vanleer,dulo,du,h);
      }
      // Compute error. LOfu - L^{exact}(u)_i
      LOfu -= LOfuE;
      errorMax[lev] = LOfu.absMax(0);
      //errorL1[lev]= norm(LOfu,h,1.0,0);
      errorL1[lev] = LOfu.sumAbs(0)*pow(h,DIM);
      cout << "max error = " << errorMax[lev] << endl;
      cout << "L1 error = " << errorL1[lev] << endl;
      HDF5Handler h5;
      h5.writePatch({},h,LOfu,"error"+to_string(nx)+"_"+to_string(dir));
      nx*=2;
    }
  double rateMax = log(errorMax[0]/errorMax[1])/log(2.0);
  cout << "max rate (should be 1) = " << rateMax << endl;
   double rateL1 = log(errorL1[0]/errorL1[1])/log(2.0);
  cout << "L1 rate (should be 2) = " << rateL1 << endl;
   PR_TIMER_REPORT();
}                     
