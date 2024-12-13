#include "Proto.H"
#define NUMCOMPS 1
#define NGHOST 4
#define MAXLEV 2
#define ORDER 4
#define REFRATIO 4
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
  PR_TIMER_SETFILE("FVInterpolation" + to_string(nx)+".time.table");
  PR_TIMERS("main");
  Array<double,MAXLEV> error;
  for (int lev = 0; lev < MAXLEV;lev++)
    {
      // initialize.
      Box bx(Point::Ones(-NGHOST),Point::Ones(nx+NGHOST));
      Box bx0(Point::Zeros(),Point::Ones(nx-1));
      Box bx0fine = bx0.refine(Point::Ones(REFRATIO));
      
      BoxData<double,NUMCOMPS> upoint(bx);
      BoxData<double,NUMCOMPS> ufinepoint(bx0fine.grow(Point::Ones(1)));
      
      double h = 1.0/nx;
      Array<double,DIM> offset = {.5,.5,.5};
      forallInPlace_p(initialize,upoint,offset,h);
      double hfine = h/REFRATIO;
      forallInPlace_p(initialize,ufinepoint,offset,hfine);
      BoxData<double,NUMCOMPS> ufine = Operator::convolve(ufinepoint);
      cout << ufine.box() << endl;
      BoxData<double,NUMCOMPS> u = Operator::convolve(upoint);
      
      // define Interpolation operator.
      auto I = InterpStencil<double>::FiniteVolume(Point::Ones(REFRATIO),5);
      //auto I = InterpStencil<double>::Constant(Point::Ones(4));
      
      // apply operator.      {
      BoxData<double,NUMCOMPS> IOfu(bx0fine);
        {
          PR_TIMERS("timed section");
          IOfu |= I(u,1.0);
        }
      // Compute error I(u) - ufine.
        HDF5Handler h5;
        BoxData errbd(ufine.box());
        IOfu.copyTo(errbd);
        cout << errbd.box() << endl;
        h5.writePatch({},hfine,IOfu,"IOfu"+to_string(nx));
        h5.writePatch({},hfine,ufine,"uexact"+to_string(nx));
        errbd *= -1;
        errbd += ufine;
        error[lev] = errbd.absMax(0);
        cout << "max error  = " << error[lev] << endl;
        h5.writePatch({},hfine,errbd,"error"+to_string(nx));
       
        nx*=2;
    }
  double rate = log(error[0]/error[1])/log(2.0);
  cout << "rate (should be 4) = " << rate << endl;
   PR_TIMER_REPORT();
}                      
