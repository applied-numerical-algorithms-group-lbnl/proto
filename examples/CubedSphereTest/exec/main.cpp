#include "Proto.H"
#include "TestMapping.H"

void GetCmdLineArgumenti(int argc, const char** argv, const char* name, int* rtn)
{
    size_t len = strlen(name);
    for(int i=1; i<argc; i+=2)
    {
        if(strcmp(argv[i]+1,name) ==0)
        {
            *rtn = atoi(argv[i+1]);
            std::cout<<name<<"="<<" "<<*rtn<<std::endl;
            break;
        }
    }
}
int main(int argc, char* argv[])
{
    // Setting up simple example of Poisson solver where the
    // input and output data is allocated by the user.
  PROTO_ASSERT(DIM==3,"cubed-sphere works only for DIM=3");
  int nx = 8;
    GetCmdLineArgumenti(argc, (const char**)argv, "nx", &nx);
    std::cout << "command line input:" << std::endl;
    std::cout << argv[0] << " -nx " << nx  << std::endl << endl;
    PR_TIMER_SETFILE(to_string(nx) + ".DIM" + to_string(DIM) + ".cubedSphere.time.table");
    PR_TIMERS("main");
    int nGhost = 4;  
    int maxRef = 1;
    Array<BoxData<double>,3> JRefs;
    HDF5Handler h5;
    for (int ref=0;ref < maxRef; ref++)
      {
        BoxData<double,DIM> X;
        Array<BoxData<double,DIM>,DIM> NT;
        Box bx(Point::Zeros(),Point::Ones(nx-1));
        BoxData<double>& J = JRefs[ref];
        Box bxnode = bx.extrude(Point::Ones());
        BoxData<double> radius(bxnode);
        double h = 1.0/nx;
        forallInPlace_p([]PROTO_LAMBDA(Point a_pt,                            
                                       Var<double,1>& a_rad,
                                       double a_h)
                        {
                          double xi0 = 1.0 + a_pt[0]*a_h;
                          a_rad(0) = xi0*xi0/2;
                        },bxnode,radius,h);
        /*
          Box bx0 = bx.extrude(Point::Basis(0));
          BoxData<double,DIM> nsphere(bxnode);
          Box bx02 = bx.extrude(Point::Ones() - Point:Basis(1));
          BoxData<double,DIM> n_x_d2n(bx01);  
          Box bx01 = bx.extrude(Point::Ones() - Point:Basis(2));
          BoxData<double,DIM> d1n_x_n(bx02);
          cout << bx0 << endl;
          cout << bx02 << endl;
          cout << bx01 << endl;
          cubedSphere_nSphere(nsphere,h);
          cubedSphere_edge1(d1n_x_n,h);
          cubedSphere_edge2(n_x_d2n,h);
          BoxData<double,DIM> nsphereFaceA = Stencil<double>::cornersToFaces(0)(nsphere);
          BoxData<double,DIM> d2_d1n_x_n = Stencil<double>::Fluxdivergence(2)(d1n_x_n,-.5);
        */
        Operator::cubedSphereGeometry(X,NT,J,radius,bx,h);
        h5.writePatch(h,X,"map"+to_string(nx));
        for (int dir = 0; dir < DIM; dir++)
          {
            NT[dir]*=nx*nx;
            h5.writePatch(h,NT[dir],"cofactor_"+to_string(dir)+"_"+to_string(nx));
          }
        J *= nx*nx;
        BoxData<double,DIM> divNT(bx);
        divNT.setToZero();
        for (int dir = 0; dir < DIM; dir++)
          {
            divNT += Stencil<double>::FluxDivergence(dir)(NT[dir]);
          }
        cout << endl << "divNT:"<<endl;
        divNT.print();
        h5.writePatch(h,divNT,"divNT"+to_string(nx));
        cout << endl << "Jacobian:"<<endl;
        h5.writePatch(h,JRefs[ref],"jacobian"+to_string(nx));
        JRefs[ref].print();
        nx*=2;
      }
    
    if (maxRef >= 3)
      {
        BoxData<double> JAv;
        JAv = Stencil<double>::AvgDown(2)(JRefs[1],-1.0);
        JAv+=JRefs[0];
        cout <<endl<<"Jerr0 norm: "<< JAv.absMax() <<endl;
        h5.writePatch(4./nx,JAv,"Jerr0"+to_string(nx));
        JAv.print();
        JAv = Stencil<double>::AvgDown(2)(JRefs[2],-1.0);
        JAv+=JRefs[1];
        h5.writePatch(2./nx,JAv,"Jerr1"+to_string(nx));
        cout <<endl<<"Jerr1 norm: "<< JAv.absMax() <<endl;
        JAv.print();
      }
    
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

