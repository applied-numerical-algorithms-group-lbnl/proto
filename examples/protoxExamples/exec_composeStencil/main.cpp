#include "Proto.H"
#define NUMCOMPS 1
#define NGHOST 1
#define MAXLEV 2
#define ORDER 4
using namespace Proto;
template<typename T,unsigned int C, MemType MEM>
PROTO_KERNEL_START
void initialize_F(Point a_pt,
                  Var<T,C,MEM>& a_u,
                  Array<T,DIM> a_offset,
                  T a_h)
{
  for (int comp = 0; comp < NUMCOMPS;comp++)
    {
      a_u(comp) = 1.;
      for (int dir = 0; dir < DIM; dir++)
        {
          a_u(comp) *= sin(a_pt[dir]*M_PI*a_h + a_offset[dir]*M_PI*a_h);
        }
    }
}
PROTO_KERNEL_END(initialize_F,initialize)
template<typename T,unsigned int C, MemType MEM>
T norm(
       BoxData<T,C,MEM>& a_phi,
       T a_h,
       T a_order,
       int a_comp)
{
  BoxData<T,1,MEM> absuToTheP =
    forall<T,1,MEM>([ ] PROTO_LAMBDA(Var<T,1,MEM>& a_output,
                            const Var<T,C,MEM>& a_phi,
                            T a_h,
                            T a_order,
                            int a_comp)
           {
             a_output(0) = pow(abs(a_phi(a_comp)),a_order)
             *pow(a_h,DIM);
           },a_phi,a_h,a_order,a_comp);
  T norm = pow(absuToTheP.sum(0),1.0/a_order);
  return norm;
}
int main(int argc, char *argv[])
{
  
  int nx;
  cout << "input nx" << endl;
  cin >> nx;
  PR_TIMER_SETFILE("composeStencils" + to_string(nx)+".time.table");
  PR_TIMERS("Laplacian_stencil");
  Array<double,MAXLEV> error;
  for (int lev = 0; lev < MAXLEV;lev++)
    {
      // initialize.
      Box bxcorner0(Point::Zeros(),Point::Ones(nx));
      Box bxcenter0(Point::Zeros(),Point::Ones(nx-1));
      Box bxcorner = bxcorner0.grow(Point::Ones(NGHOST));
      Array<double,DIM> offsetCorner = {0.,0.,0.};
      Array<double,DIM> offsetCenter = {.5,.5,.5};
      double h = 1.0/nx;
      
      BoxData<double,NUMCOMPS> ucorner =
        forall_p<double,NUMCOMPS>(initialize,bxcorner,offsetCorner,h);
      BoxData<double,NUMCOMPS> ucenter =
        forall_p<double,NUMCOMPS>(initialize,bxcenter0,offsetCenter,h);
      
      // define stencils.
      Stencil<double> I2 = 1.0*Shift(Point::Zeros()) + 1.0*Shift(Point::Basis(0,1));
      for (int dir = 1; dir < DIM; dir++)
        {
          I2 = I2 * (1.0*Shift(Point::Zeros()) + 1.0*Shift(Point::Basis(dir,1)));
        }
      I2 = (1.0/8.0)*I2;
      Stencil corr = (1.0*Shift(Point::Zeros())
                      + (-1.0/8.0)*Stencil<double>::Laplacian());
 
      // apply operator.
      BoxData<double,NUMCOMPS> IOfu;
      {
        PR_TIMERS("timed section");
        BoxData<double,NUMCOMPS> u1 = I2(ucorner);
        IOfu = corr(u1);
      }
   
      // Compute error. 
      BoxData<double,NUMCOMPS> errbd(bxcenter0);
      ucenter.copyTo(errbd);
      errbd *=-1;
      errbd += IOfu;
      error[lev] = errbd.absMax(0);
      cout << "max error = " << error[lev] << endl;
      HDF5Handler h5;
      h5.writePatch({},h,errbd,"error"+to_string(nx));
      nx*=2;
    }
  double rate = log(error[0]/error[1])/log(2.0);
  cout << "rate = " << rate << endl;
   PR_TIMER_REPORT();
}
         
                        
