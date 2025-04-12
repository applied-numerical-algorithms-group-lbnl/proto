#include <iostream>
#include <iomanip>
#include "Proto.H"
#include "Multigrid.H"
#include "InputParser.H"
#include "MakeLaplaceStencils.H"

using namespace Proto;

int TIME_STEP = 0;
int PLOT_NUM = 0;
PROTO_KERNEL_START
void
f_initializeF(
              const Point& a_pt,
              Var<double>& a_rhs,
              const double& a_h)
{
        a_rhs(0) = 0.;
}
PROTO_KERNEL_END(f_initializeF, f_initialize);

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    using Proto::pout; 

    int domainSize = 128;
    int boxSize = 32;
    int numLevels = 7;
    double tolerance = 1.e-15;
    int solverType = 0;
    int maxIter = 20;
    
    double t0 = 0.;
    int maxRefs = 2;
    //std::array<bool, DIM> periodicity;
    Array<bool, DIM> periodicity;
    periodicity.fill(false);

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("boxSize",         boxSize);
    args.add("numLevels",       numLevels);
    args.add("maxIter",         maxIter);
    args.add("tolerance",       tolerance);
    args.add("solverType",      solverType);
       
    args.parse(argc, argv);
    if (procID() == numProc()-1)
      {
        cout << procID() << endl;
        args.print();
      }
    Point boxSizeVect = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);

    double physDomainSize = 1.0;
    PR_TIMER_SETFILE("Multigrid"  + to_string(DIM) + "D" +"_" + to_string(boxSize) + "_" + to_string(domainSize)+".time.table");
    Array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);
    Multigrid<double,HOST> mg(layout,solverType,numLevels);
  
    Stencil<double> sten = LaplaceStencil<double>(solverType);
    Stencil<double> lap2 = Stencil<double>::Laplacian();
    int order = 2;

    if (solverType > 0) order = 4;
    int diamStencil = sten.ghost()[0];
    LevelBoxData<double,1,HOST> phi(layout,Point::Ones(diamStencil));
    LevelBoxData<double,1,HOST> rhs(layout,Point::Zeros());
    LevelBoxData<double,1,HOST> rhs0(layout,Point::Zeros());

    for (auto dit : layout)
      {
        forallInPlace_p(f_initialize,layout[dit],rhs0[dit],dx[0]);
        phi[dit].setVal(0.);
      }
    Reduction<double,Abs,HOST> rxn0;
    for (auto dit : layout)
      {
        double maxres = rhs[dit].absMax();
        rxn0.reduce(&maxres,1);
      }
    double resnorm0 = rxn0.fetch();
    if (procID() == 0) cout << "Initial residual = " << resnorm0 << endl;
    h5.writeLevel({},dx,rhs0,"rhsInitial" + to_string(DIM) + "D");
    h5.writeLevel({},dx,rhs,"resInitial" + to_string(DIM) + "D");
    h5.writeLevel({},dx,phi,"phiInitial" + to_string(DIM) + "D");
    // Begin multigrid iteration.
   
    for (int iter = 0; iter < maxIter; iter++)
      {
        // FAS input.
        mg.residual(rhs,phi,rhs0,dx[0]);
        for (auto dit : layout)
          {
            rhs[dit] += mg.smoother()(phi[dit],1.0/(dx[0]*dx[0]));
          }
        
        mg.mgRelax(phi,rhs0,dx[0],numLevels);
        mg.residual(rhs,phi,rhs0,dx[0]);
        
        h5.writeLevel({},dx,phi,"phi" + to_string(DIM) + "D_"+to_string(iter));
        h5.writeLevel({},dx,rhs,"res" + to_string(DIM) + "D_"+to_string(iter));
       
        Reduction<double,Abs,HOST> rxn;
          for (auto dit : layout)
             {
               double maxres = rhs[dit].absMax();
               rxn.reduce(&maxres,1);
             }
          double resnormIter = rxn.fetch();
        if (procID() == 0) cout << "residual at iteration "
                                << iter << " = " << resnormIter << endl;
        if (resnormIter < resnorm0*tolerance) break;
      }
    h5.writeLevel({},dx,rhs,"resFinal" + to_string(DIM) + "D");
    h5.writeLevel({},dx,phi,"phiFinal" + to_string(DIM) + "D"); 
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
