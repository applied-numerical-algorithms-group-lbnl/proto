#include "Proto.H"
#include "Proto_LevelBoxData.H"
#include "Proto_DisjointBoxLayout.H"
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

PROTO_KERNEL_START void fusedRelaxT(Var<double> u, Var<double> temp, Var<double> rhs, double lambda)
{
  u(0) += lambda*(rhs(0)-temp(0));
}
PROTO_KERNEL_END(fusedRelaxT, fusedRelax);

#define PFLOPS 3
  
PROTO_KERNEL_START void initRHST(const Point& p, Var<double> rhs)
{
  rhs(0) = 0.0;
  if(p==Point(8,8,8)) rhs(0) = 1.0;
  if(p==Point(12,12,8)) rhs(0) = -1.0;
}

PROTO_KERNEL_END(initRHST, initRHS);

void InitializeMMBData( 
                       MMBLevelBoxData<double>& a_phi,
                       const MMBMaps& a_mmbmap
                        )
{
  Stencil<double> cToF;
  for (int dir = 0;dir < DIM;dir++)
    {
      cToF *= .5*Shift::Basis(dir,1) + .5*Shift::Zeros();
    }
  cToF += (1./12.)*Stencil<double>::Laplacian();
  for (auto dit = a_phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double,DIM> X = a_mmbmap.map(*dit,Point::Ones(4));
      BoxData<double,DIM> Xcen = cToF(X);
      auto initPointData = forall_p(initFunction,Xcen);
      a_phi[*dit] = Stencil<double>::convolution(4)(initPointData);
    }
}
int main(int argc, char* argv[])
{
  
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");

  // Create MMBPDs
  MMBProblemDomain mmpbd;
  MMBMaps mmbmaps;
  InitializeMMB(mmbpd,mmbmaps); // provided by the application.
  vector<int> blocksizes = mmbpd.defaultBlocksizes();
  // Create MMBDBLs.  
  MMBDisjointBoxLayout mmbdbl(mmbpd,blocksizes);

  /* Geometry data holders for mapping. We will compute once and store 
     the Jacobian at cell centers (J), cofactor matrices (NT) on faces. */

  // Create and fill containers for J, normal components of cofactor matrices.
  
  MMBLayoutData<array<BoxLayout<double,DIM>,DIM> > NT(mmbdbl);
  MMBLayoutData<BoxData<double> > J(mmbdbl);
  for (auto dit = mmdbl.begin();*dit != mmdbl.end();++dit)
    {
      auto &Xbox = X[*dit];
      auto &NTBox = NT[*dit];
      for (int dir = 0; dir < DIM; dir++)
        {
          NTBox[dir] = cofactor(X,dir);
        }
      auto &Jbox = J[*dit];
      JBox = jacobian(X,NT);
    }
     
  // Create MMB LevelBoxData.
  array<Point,4> ghostSizes;
  ghostsizes.fill(Point::Ones(3));
  
  // ghostsizes[0] = interior, ghostsizes[d] = codimension d+1 block boundary.
  MMBLevelBoxData phi(mmbdbl,ghostSizes);

  // Initialize phi.
  MMBInitialize(phi,mmbmap);

  // Evaluate Laplacian(phi) = div(NT.((1/J)*N).grad(phi)).
  
  for (auto dit = mmdbl.begin();*dit != mmdbl.end();++dit)
    {
      // Calculation of 
      auto &Xbox = X[*dit];
      auto &Jbox = J[*dit];
      BoxData<double>& phibox = phi.mmbInterpolate(*dit,5);
      for (int dir = 0; dir < DIM; dir++)
        {
          BoxDataMatrix NTMatrix(NT[*dit][dir].box());
          NTMatrix.setToZero();
          for (int dir1 = 0; dir1 < DIM; dir++)
            {
              for (int dir2 = 0;dir2 < DIM; dir2++)
                {
                  auto NTcomp = slice(NTbox[dir2],dir1);
                  auto NTMatrixcomp = slice(NTMatrix,dir1,dir2);
                  if (dir1==dir) 
                    {
                      NTcomp.copyTo(NTMatrixComp);
                    }
                  else
                    {
                      NTMatrixComp +=
                        Stencil<double>::FaceAvToFaceAv(NTcomp,dir1,dir2,4);
                    }
                }
            }
          BoxData<double> JFace = Stencil<double>::CellToFace(dir,Side::Lo,4)(Jbox);

          
          BoxData<double,DIM> gradPhiFace4(NTBox.box());
          BoxData<double,DIM> gradPhiFace2(NTBox.box());

          // Compute (1/J)*NT.
          auto JinvNTFace = faceMatrixQuotient(J,NTFace,J,NTFace,dir);
          // compute gradient on face.
          for (int dir2 = 0; dir < DIM; dir++)
            {
              auto &gradPhiComp4 = slice(gradPhiFace4,dir,1,1);
              auto &gradPhiComp2 = slice(gradPhiFace2,dir,1,1);
              if (dir2 == dir)
                {
                  gradPhiComp4 = Stencil<double>::diffCellToFace(dir2)(phi);
                  gradPhiComp2 = Stencil<double>::diffCellToFace(dir2,2)(phi);
                }
              else
                {
                  gradPhiComp4 =
                    (Stencil<double>::diffCellToFace(dir2)
                    *Stencil<double>::cellToFace(dir))(phi);
                  gradPhiComp2 =
                    (Stencil<double>::diffCellToFace(dir2,2)
                     *Stencil<double>::cellToFace(dir,2))(phi);
                }              
            }
          // Product of transpose(NT) and the gradient.
          auto gradxPhiFace4 = faceMatrixProductAB
            (gradphiFace4,JinvNTBox,gradFace2,JinvNTBox,dir);
          auto gradPhiFace2 = matrixProduct2(gradPhi2,JinvNTBox);
          // dir component of NT.gradxPhi.
          auto flux = faceMatrixProductATB
            (NTBox,gradxPhi4,NTBox,gradxPhi);
          LaplacianPhiBox += Stencil<double>::divergence(dir)(flux);
          LaplacianBox /= JBox;
        }      
  PR_TIMER_REPORT();
  return 0;
}
