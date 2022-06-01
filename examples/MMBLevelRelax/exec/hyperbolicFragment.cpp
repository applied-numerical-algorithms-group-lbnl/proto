MMBProblemDomain mmpbd;
MMBMaps mmbmaps;
InitializeMMB(mmbpd,mmbmaps); // provided by the application.
vector<Point> blocksizes = mmbpd.defaultBlocksizes();
  
MMBDisjointBoxLayout mmbdbl(mmbpd,blocksizes);

// Create and fill containers for J, normal components of cofactor matrices.
  
MMBLayoutData<array<BoxData<double,DIM>,DIM> > NT(mmbdbl);
MMBLayoutData<BoxData<double> > J(mmbdbl);
for (auto dit = mmdbl.begin();*dit != mmdbl.end();++dit)
  {
    BoxData<double,DIM> X = mmbmaps.map(*dit,Point::Ones(4));
    auto &NTbox = NT[*dit];
    for (int dir = 0; dir < DIM; dir++)
      {
          NTbox[dir] = cofactor(X,dir);
      }
    auto &Jbox = J[*dit];
    JBox = jacobian(Xbox,NTbox);
  }
// Create and initialize MMB LevelBoxData.

array<Point,DIM+1> ghostSizes;
ghostsizes.fill(Point::Ones(3));

// ghostsizes[0] = interior, ghostsizes[d] = codimension d+1 block boundary.

MMBLevelBoxData<double> phi(mmbdbl,ghostsizes);

for (auto dit = phi.begin();*dit != dit.end();++dit)
    {
      BoxData<double,DIM> X = a_mmbmap.map(*dit,Point::Ones(4));
      BoxData<double,DIM> Xcen = Stencil<double>::cornerToCenter()(X);
      auto initPointData = forall_p(initFunction,Xcen);
      phi[*dit] = Stencil<double>::convolution(4)(initPointData);
    }

// Compute the divergence of the fluxes.

phi.exchange();
MMBLevelBoxData<double> divF(mmbdbl);
for (auto dit = phi.begin();*dit != dit.end();++dit)
    {
     BoxData<double>& phibox = phi.mmbInterpolate(*dit,4);
     BoxData<double>& divFBox = divF[*dit];
     for (int dir = 0; dir < DIM; dir++)
        {
          //...
          auto flux =
            Operator<double>::faceMatrixProductATB(NTBox,FluxAv4,NTBox,FluxAv2,dir);
          divFBox += Stencil<double>::divergence(dir)(flux);
        }
    }
// Missing flux matching steps.
