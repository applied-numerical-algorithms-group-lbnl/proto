#include "Proto.H"
#define NUMCOMPS 1
#define NGHOST 2
#define MAXLEV 2
#define ORDER 4
using namespace Proto;
int main(int argc, char *argv[])
{
  
  int nx;
  cout << "input nx" << endl;
  cin >> nx;
  PR_TIMER_SETFILE("laplacian_stencil" + to_string(nx)+".time.table");
  PR_TIMERS("Laplacian_stencil");
  Array<double,MAXLEV> error;
  for (int lev = 0; lev < MAXLEV;lev++)
    {
      // initialize.
      Box bx(Point::Ones(-NGHOST),Point::Ones(nx+NGHOST));
      Box bx0(Point::Ones(1),Point::Ones(nx-1));
      BoxData<double,NUMCOMPS> u(bx);
      double h = 1.0/nx;
      forallInPlace_p([ ]PROTO_LAMBDA(Point a_pt,
                                      Var<double,NUMCOMPS>& a_u,
                                      double a_h)
                      {
                        for (int comp = 0; comp < NUMCOMPS;comp++)
                          {
                            a_u(comp) = 1.;
                            for (int dir = 0; dir < DIM; dir++)
                              {
                                a_u(comp) *= sin(a_pt[dir]*M_PI*a_h);
                              }
                          }
                      },bx,u,h);
      
      // define operator.
      Stencil<double> L;
      for (int dir = 0; dir < DIM; dir++)
        {
          L += Stencil<double>::Derivative(2,dir,ORDER);
        }      
      // apply operator.
      BoxData<double,NUMCOMPS> LOfu;
      {
        PR_TIMERS("timed section");
        LOfu = L(u,bx0,1/(h*h));
      }
      // Compute error. L^{exact}(u)_i = -pi*pi*DIM*u_i
      BoxData<double,NUMCOMPS> LOfuE(u.box());
      u.copyTo(LOfuE);
      LOfuE *= (-M_PI*M_PI*DIM);
      LOfu -= LOfuE;
      error[lev] = LOfu.absMax(0);
      cout << "max error = " << error[lev] << endl;
      HDF5Handler h5;
      h5.writePatch({},h,LOfu,"error"+to_string(nx));
      nx*=2;
    }
  double rate = log(error[0]/error[1])/log(2.0);
  cout << "rate = " << rate << endl;
   PR_TIMER_REPORT();
}
         
                        
