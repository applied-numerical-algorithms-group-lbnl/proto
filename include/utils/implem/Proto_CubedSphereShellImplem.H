MBProblemDomain CubedSphereShell::Domain(
    int a_domainSize,
    int a_thickness,
    int a_radialDir)
{
    MBProblemDomain domain(6);
    auto CCW = CoordPermutation::ccw(a_radialDir);
    auto CW = CoordPermutation::cw(a_radialDir);
    CoordPermutation RNorth, RSouth, RI;
    Point x = Point::X();
    Point y = Point::Y();
    Point z = Point::Z();
    for (int bi = 2; bi < 6; bi++)
    {
        int srcBlock = bi;
        int dstBlock = bi+1;
        if (dstBlock > 5) { dstBlock = 2; }
        switch (a_radialDir)
        {
            case 0:
                domain.defineBoundary(srcBlock, dstBlock, y, RI);
                domain.defineBoundary(srcBlock, 1, z, RNorth);
                domain.defineBoundary(srcBlock, 0, -z, RSouth);
                break;
            case 2:
                domain.defineBoundary(srcBlock, dstBlock, x, RI);
                domain.defineBoundary(srcBlock, 1, y, RNorth);
                domain.defineBoundary(srcBlock, 0, -y, RSouth);
                break;
            default:
                MayDay<void>::Abort("NOT IMPLEMENTED");
        }
        RNorth = RNorth*CW;
        RSouth = RSouth*CCW;
    }
    domain.graph().fixRotations();
    Point domainSize = Point::Ones(a_domainSize);
    domainSize[a_radialDir] = a_thickness;
    for (int bi = 0; bi < 6; bi++)
    {
        domain.defineDomain(bi, domainSize);
    }
    return domain;
}
    
template<MemType MEM>
MBLevelMap<MBMap_CubedSphereShell, MEM>
CubedSphereShell::Map(
        const MBDisjointBoxLayout&  a_layout,
        Point   a_ghost)
{
    MBLevelMap<MBMap_CubedSphereShell, HOST> map;
    map.define(a_layout, a_ghost);
    return map;
}

template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBLevelMap<MBMap_CubedSphereShell, MEM>
CubedSphereShell::Map(
        const MBLevelBoxData<T,C,MEM,CTR>& a_data)
{
    MBLevelMap<MBMap_CubedSphereShell, HOST> map;
    map.define(a_data.layout(), a_data.ghost());
    return map;
}

template<template<typename, template<MemType> typename, MemType> class OpType,
    typename T,
    MemType MEM>
MBLevelOp<OpType, MBMap_CubedSphereShell, T>
CubedSphereShell::Operator(
    MBLevelMap<MBMap_CubedSphereShell, MEM>& a_map)
{
    MBLevelOp<OpType, MBMap_CubedSphereShell, T> op;
    op.define(a_map);
    return op;
}

template<MemType MEM>
MBInterpOp
CubedSphereShell::InterpOp(
        const MBDisjointBoxLayout& a_layout,
        Point a_ghost,
        int a_order)
{
    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(2))
    {
        if (pi.abs().sum() <= 2 
            && pi[CUBED_SPHERE_SHELL_RADIAL_COORD] == 0)
        {
            footprint.push_back(pi);
        }
    }
    MBInterpOp op(a_ghost, a_order);
    op.copyAxis(CUBED_SPHERE_SHELL_RADIAL_COORD, true);
    for (int bi = 0; bi < 6; bi++)
    {
        MBLevelMap<MBMap_CubedSphereShellPolar, MEM> polarMap;
        polarMap.define(a_layout, a_ghost, bi);
        op.define(polarMap, footprint, a_order, bi); 
    }
    return op; 
}
    
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBInterpOp
CubedSphereShell::InterpOp(
        MBLevelBoxData<T,C,MEM,CTR>& a_data,
        int a_order)
{
    return InterpOp<MEM>(a_data.layout(), a_data.ghost()[0], a_order);
}

template<typename T>
int _findIndex(std::vector<T> sorted, T value)
{
    for (int ii = 0; ii < sorted.size(); ii++)
    {
        if (sorted[ii] > value) { return ii; }
    }
    return sorted.size()-1;
}

template<typename T>
Array<T, 3> _getPolar(Array<T,3> X)
{
    T x = X[0]; T y = X[1]; T z = X[2];
    int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
    int thetaCoord = (rCoord + 1) % 3;
    int phiCoord = (rCoord + 2) % 3;
    
    Array<T,3> R;
    T r = sqrt(x*x + y*y + z*z);
    R[rCoord] = r;
    R[thetaCoord] = acos(z/r);
    T phi = atan2(y,x);
    if (phi < 0) { phi += 2*M_PI; }
    R[phiCoord] = phi;
    return R;
}
#if 0
Point _getPolarImage(Box b, Point p)
{

    //theta = coordinate 0;
    //phi = coordinate 1;
    //origin at north pole
    Point shift = b.low();
    Box b0 = b.shift(-shift);
    Point q = p - shift;
    
    int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
    int o = (rCoord + 1) % 3;
    int l = (rCoord + 2) % 3;
    
    while (q[l] < 0)
    {
        q[l] += b0.sizes()[l];
    }
    q[l] = q[l] % b0.sizes()[l];
    if (q[o] < 0)
    {
        q[l] = (q[l] + b0.sizes()[l]/2) % b0.sizes()[l];
        q[o] = -(1 + q[o]);
    } else if (q[o] >= b0.sizes()[o])
    {
        q[l] = (q[l] + b0.sizes()[l]/2) % b0.sizes()[l];
        q[o] = b0.sizes()[l] - (q[o] % b0.sizes()[l]) - 1;
    }
    return q + shift;
}
#endif

template<typename T, MemType MEM>
MBInterpOp
CubedSphereShell::BCInterpOp(
        const MBLevelMap<MBMap_CubedSphereShell, MEM>& a_map,
        MBDisjointBoxLayout&       a_srcLayout,
        std::vector<T>             a_dtheta,
        Side::LoHiSide             a_side)
{
    HDF5Handler h5;
    // PROTO_ASSERT(a_srcLayout.numBlocks() == 1,
    //         "CubedSphereShell::BCInterpOp | Input data should be on a single block layout.");

    // SETUP
    const auto& layout = a_map.layout();
    int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
    int thetaCoord = (rCoord + 1) % 3;
    int phiCoord = (rCoord + 2) % 3;
    std::vector<string> varnames(3);
    varnames[rCoord] = "R";
    varnames[thetaCoord] = "Theta";
    varnames[phiCoord] = "Phi";

    int rSize = layout.domain().getBlock(0).box().sizes()[rCoord];
    Point srcBoxSizes;
    srcBoxSizes[rCoord] = 1;
    srcBoxSizes[thetaCoord] = a_srcLayout.domain().getBlock(0).box().sizes()[thetaCoord];
    srcBoxSizes[phiCoord] = a_srcLayout.domain().getBlock(0).box().sizes()[phiCoord];
    Box srcBox(srcBoxSizes);
    T R = 1.0;
    int RIndex = 0;
    if (a_side == Side::Lo)
    {
        R = CUBED_SPHERE_SHELL_R0; 
        srcBox = srcBox.shift(Point::Basis(rCoord, -1));
        RIndex = -1;
    } else {
        R = CUBED_SPHERE_SHELL_R1; 
        srcBox = srcBox.shift(Point::Basis(rCoord, rSize));
        RIndex = rSize;
    }
    
    // CACHE THETA ON EACH LATITUDE (used for fast search)
    std::vector<T> thetaArray(srcBoxSizes[thetaCoord]);
    for (int ii = 0; ii < srcBoxSizes[thetaCoord]; ii++)
    {
        thetaArray[ii] = (ii == 0 ? 0 : thetaArray[ii-1] + a_dtheta[ii-1]);
    }

    // BUILD POLAR MAP ON SRC BOX
    BoxData<T,3,MEM> srcMap(srcBox);
    T dphi = 2.0*M_PI / srcBoxSizes[phiCoord];
    for (auto pi : srcBox)
    {
        srcMap(pi,rCoord) = R;
        srcMap(pi,thetaCoord) = thetaArray[pi[thetaCoord]];
        srcMap(pi,phiCoord) = (pi[phiCoord] + 0.5)*dphi;
    }
    h5.writePatch(varnames, srcMap, "SRC_MAP");


    // BUILD THE INTERPOLATION OPERATOR
    MBInterpOp op;
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        // All points in this block (not including ghost)
        Box blockDomain = layout.domain().getBlock(block).box();
        // Boundary of this block and this patch respectively
        Box boundDomain, boundDomain_i;

        if (a_side == Side::Lo)
        {
            boundDomain = blockDomain.face(rCoord, a_side, 1);
            boundDomain_i = layout[iter].face(rCoord, a_side, 1) & boundDomain;
        } else {
            boundDomain = blockDomain.adjacent(rCoord, a_side, 1);
            boundDomain_i = layout[iter].adjacent(rCoord, a_side, 1) & boundDomain;
        }
        if (boundDomain_i.size() <= 0) { continue; }

        // Get Cartesian coordinates of cubed sphere map (J is unused; ignore it)
        BoxData<double, DIM, MEM> X(boundDomain_i);
        BoxData<double, 1, MEM> J(boundDomain_i);
        a_map.apply(X, J, block);
        if (a_side == Side::Lo)
        {
            X.shift(Point::Basis(rCoord, Side::Lo));
        }
        // CREATE POINT INTERP OP FOR EACH POINT IN CUBE SPHERE BOUNDARY
        for (auto p_dst : X.box())
        {
            auto x_dst = X.array(p_dst);    //get cartesian coords at this point
            auto r_dst = _getPolar(x_dst);  //convert cartesian to spherical polar
            
            // FIND POINT P_SRC IN POLAR DATA CLOSE TO P_DST IN ANGULAR DISTANCE
            Point p_src = Point::Zeros();
            
            p_src[rCoord] = RIndex;
            p_src[phiCoord] = r_dst[phiCoord] / dphi; //works because dphi is constant
            p_src[thetaCoord] = _findIndex(thetaArray, r_dst[thetaCoord]); // same idea but accounts for non constant dtheta
            auto r_src = srcMap.array(p_src);

            // FIND DIRECTION FROM P_SRC TO P_DST
            // NB: dr can be very close to zero, so don't divide by it
            auto dr = r_dst - r_src;
            Point s0 = dr[thetaCoord] < 0 ? Point::Basis(thetaCoord,-1) : Point::Basis(thetaCoord);
            Point s1 = dr[phiCoord] < 0 ? Point::Basis(phiCoord,-1) : Point::Basis(phiCoord);
            // leaving srcBox is forbidden due to polar singularity
            if (!srcBox.contains(p_src + s0)) { s0 = -s0; }
            if (!srcBox.contains(p_src + s1)) { s1 = -s1; }
            
            // SHIFTS DEFINE 4 POINTS FOR BILINEAR INTERPOLATION
            std::vector<Point> srcPoints;
            srcPoints.push_back(p_src);
            srcPoints.push_back(p_src + s0);
            srcPoints.push_back(p_src + s1);
            srcPoints.push_back(p_src + s0 + s1);

            // GET COORDINATES OF INTERPOLATING POINTS
            std::vector<Array<T,3>> srcData(4);
            for (int ii = 0; ii < 4; ii++)
            {
                srcData[ii] = srcMap.array(srcPoints[ii]);
            }

            // COMPUTE COEFFICIENTS
            Array<T,2> y{r_dst[thetaCoord], r_dst[phiCoord]};
            Array<T,2> y0{srcData[0][thetaCoord], srcData[0][phiCoord]};
            Array<T,2> y1{srcData[3][thetaCoord], srcData[3][phiCoord]};
            auto dy = y1 - y0;            

#if 1
            Matrix<T> A1({{1, y[0], y[1], y[0]*y[1]}});
            Matrix<T> A2({{y1[0]*y1[1], -y1[0]*y0[1], -y1[1]*y0[0], y0[0]*y0[1]},
                         {     -y1[1],        y0[1],        y1[1],      -y0[1]},
                         {     -y1[0],        y1[0],        y0[0],      -y0[0]},
                         {          1,           -1,           -1,           1}});
            auto S = A1*A2;
            S *= (1.0/(dy[0]*dy[1]));
#else
            Matrix<T> B({{1, y[0], y[1], y[0]*y[1]}});
            Matrix<T> A({{1, y0[0], y0[1], y0[0]*y0[1]},
                         {1, y1[0], y0[1], y1[0]*y0[1]},
                         {1, y0[0], y1[1], y0[0]*y1[1]},
                         {1, y1[0], y1[1], y1[0]*y1[1]}});
            auto Ainv = A.inverse();
            auto S = B*Ainv;
#endif
            
            PROTO_ASSERT(abs(S.sum()-1) < 1e-10,
                    "CubedSphereShell::BCInterpOp | Error: Invalide Interpolation Coefficients");

            // BUILD THE OPERATOR
            std::vector<T> coefs(4);
            for (int ii = 0; ii < 4; ii++) { coefs[ii] = S(0,ii); };
            std::vector<MBDataPoint> srcs(4);
            for (auto& pi : srcPoints) { pi[rCoord] = 0; }
            MBDataPoint d_dst(iter, p_dst, layout);
            srcs[0] = MBDataPoint(*a_srcLayout.begin(), srcPoints[0], a_srcLayout);
            srcs[1] = MBDataPoint(*a_srcLayout.begin(), srcPoints[1], a_srcLayout);
            srcs[2] = MBDataPoint(*a_srcLayout.begin(), srcPoints[2], a_srcLayout);
            srcs[3] = MBDataPoint(*a_srcLayout.begin(), srcPoints[3], a_srcLayout);
           
            //arguments are dst point, src points, coefficients
            op.insert(std::make_shared<MBPointInterpOp>(d_dst, srcs, coefs));
        }
    }
    return op;
}


template<typename T, MemType MEM>
MBInterpOp
CubedSphereShell::BCNearestOp(
        const MBLevelMap<MBMap_CubedSphereShell, MEM>& a_map,
        MBDisjointBoxLayout&       a_srcLayout,
        std::vector<T>             a_dtheta,
        std::vector<T>             a_theta,
        std::vector<T>             a_phi,
        Side::LoHiSide             a_side,
        int a_ghost)
{
    HDF5Handler h5;

    // SETUP
    const auto& layout = a_map.layout();
    int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
    int thetaCoord = (rCoord + 1) % 3;
    int phiCoord = (rCoord + 2) % 3;


    int rSize = layout.domain().getBlock(0).box().sizes()[rCoord];
    Point srcBoxSizes;
    srcBoxSizes[rCoord] = 1;
    srcBoxSizes[thetaCoord] = a_srcLayout.domain().getBlock(0).box().sizes()[thetaCoord];
    srcBoxSizes[phiCoord] = a_srcLayout.domain().getBlock(0).box().sizes()[phiCoord];
    Box srcBox(srcBoxSizes);
    T R = 1.0;
    int RIndex = 0;
    if (a_side == Side::Lo)
    {
        R = CUBED_SPHERE_SHELL_R0; 
        srcBox = srcBox.shift(Point::Basis(rCoord, -1));
        RIndex = -1;
    } else {
        R = CUBED_SPHERE_SHELL_R1; 
        srcBox = srcBox.shift(Point::Basis(rCoord, rSize));
        RIndex = rSize;
    }
    
    // CACHE THETA ON EACH LATITUDE (used for fast search)
    std::vector<T> thetaArray(srcBoxSizes[thetaCoord]);
    for (int ii = 0; ii < srcBoxSizes[thetaCoord]; ii++)
    {
        thetaArray[ii] = (ii == 0 ? 0 : thetaArray[ii-1] + a_dtheta[ii-1]);
    }
    T dphi = 2.0*M_PI / srcBoxSizes[phiCoord];



    // BUILD THE INTERPOLATION OPERATOR
    MBInterpOp op;
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        // All points in this block (not including ghost)
        Box blockDomain = layout.domain().getBlock(block).box();
        // Boundary of this block and this patch respectively
        Box boundDomain, boundDomain_i;

        if (a_side == Side::Lo)
        {
            boundDomain = blockDomain.face(rCoord, a_side, 1);
            boundDomain_i = layout[iter].face(rCoord, a_side, 1) & boundDomain;
        } else {
            boundDomain = blockDomain.adjacent(rCoord, a_side, 1);
            boundDomain_i = layout[iter].adjacent(rCoord, a_side, 1) & boundDomain;
        }
        if (boundDomain_i.size() <= 0) { continue; }

        boundDomain_i = boundDomain_i.grow(a_ghost*Point::Basis(thetaCoord) + a_ghost*Point::Basis(phiCoord));
        // Get Cartesian coordinates of cubed sphere map (J is unused; ignore it)
        BoxData<double, DIM, MEM> X(boundDomain_i);
        BoxData<double, 1, MEM> J(boundDomain_i);
        a_map.apply(X, J, block);
        if (a_side == Side::Lo)
        {
            X.shift(Point::Basis(rCoord, Side::Lo));
        }
        // CREATE POINT INTERP OP FOR EACH POINT IN CUBE SPHERE BOUNDARY
        for (auto p_dst : X.box())
        {
            auto x_dst = X.array(p_dst);    //get cartesian coords at this point
            auto r_dst = _getPolar(x_dst);  //convert cartesian to spherical polar
            
            // FIND POINT P_SRC IN POLAR DATA CLOSE TO P_DST IN ANGULAR DISTANCE
            Point p_src = Point::Zeros();
            
            p_src[rCoord] = RIndex;
            p_src[phiCoord] = r_dst[phiCoord] / dphi; //works because dphi is constant
            p_src[thetaCoord] = _findIndex(thetaArray, r_dst[thetaCoord]); // same idea but accounts for non constant dtheta
            // p_src[phiCoord] = _findIndex(a_phi, r_dst[phiCoord]);; 
            // p_src[thetaCoord] = _findIndex(a_theta, r_dst[thetaCoord]); 

            std::vector<T> coefs(4);
            for (int ii = 0; ii < 4; ii++) { coefs[ii] = 0.25; };
            std::vector<MBDataPoint> srcs(4);
            std::vector<Point> srcPoints;
            srcPoints.push_back(p_src);
            for (auto& pi : srcPoints) { pi[rCoord] = 0; }
            MBDataPoint d_dst(iter, p_dst, layout);
            srcs[0] = MBDataPoint(*a_srcLayout.begin(), srcPoints[0], a_srcLayout);
            srcs[1] = MBDataPoint(*a_srcLayout.begin(), srcPoints[0], a_srcLayout);
            srcs[2] = MBDataPoint(*a_srcLayout.begin(), srcPoints[0], a_srcLayout);
            srcs[3] = MBDataPoint(*a_srcLayout.begin(), srcPoints[0], a_srcLayout);
          
            //arguments are dst point, src points, coefficients
            op.insert(std::make_shared<MBPointInterpOp>(d_dst, srcs, coefs));
        }
    }
    return op;
}


template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
CubedSphereShell::InterpBoundaries(
                                   MBLevelBoxData<T,C,MEM,CTR>& a_data,
                                   int a_order)
{
    auto op = InterpOp(a_data, a_order);
    op.apply(a_data, a_data);
}

template<typename T, unsigned int C, MemType MEM>
void CubedSphereShell::
consToSphNGEuler(                    
                 BoxData<T,C,MEM>& a_JU,
                 BoxData<T,1,MEM>& a_dVolr,
                 Box a_bx,
                 Box a_blockBox,
                 T a_dxiPerp,
                 int a_block,                 
                 unsigned int a_order)
{
  PR_TIMERS("consToSphNGEuler");
  
  //HDF5Handler h5;
  BoxData<T,C,MEM> USph(a_bx);
  USph.setVal(1.);
  // Get geometric data.
  Box bxMetric = a_bx.grow(1);
  BoxData<T,1,MEM> deta = detA<T,MEM>(bxMetric,a_dxiPerp,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);
      
  BoxData<T,DIM,MEM,DIM> deta_amat =
    CubedSphereShell::detA_A<T,MEM>(bxMetric,a_dxiPerp,a_block,0); 
  BoxData<T,DIM,MEM,DIM> invcoef(a_bx);
  invcoef.setToZero();
  forallInPlace(f_matinv3by3,invcoef,deta_amat);

  // Slice data and compute quotients.
  BoxData<T,1,MEM> JRho = slice(a_JU,0);
  BoxData<T,1,MEM> rho = slice(USph,0);      
  cellQuotientValAndInv(rho,deta,detaInv,JRho,a_bx,a_blockBox);
      
  BoxData<T,1,MEM> JRhoE = slice(a_JU,1+DIM);
  BoxData<T,1,MEM> rhoE = slice(USph,1+DIM);
  cellQuotientValAndInv(rhoE,deta,detaInv,JRhoE,a_bx,a_blockBox);
      
  BoxData<T,DIM,MEM> JUMom = slice<T,C,DIM,MEM>(a_JU,1);
  BoxData<T,DIM,MEM> USphMom = slice<T,C,DIM,MEM>(USph,1);
  cellMatrixQuotientMatAndInv(USphMom,deta_amat,invcoef,JUMom,a_bx,a_blockBox);

  BoxData<T,DIM,MEM> JUB = slice<T,C,DIM,MEM>(a_JU,2+DIM);
  BoxData<T,DIM,MEM> USphB = slice<T,C,DIM,MEM>(USph,2+DIM);
  cellMatrixQuotientMatAndInv(USphB,deta_amat,invcoef,JUB,a_bx,a_blockBox);

  // scale out radial dVol.
  BoxData<T,1,MEM> dVolrInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_dvolr)
     {
       a_output(0) = 1.0/a_dvolr(0);
     }
     ,a_dVolr); 
  for (int comp = 0; comp < C; comp++)
    {
      BoxData<T,1,MEM> USphcomp = slice(USph,comp);
      BoxData<T,1,MEM> JUcomp = slice(a_JU,comp);
      BoxData<T,1,MEM> USphtemp(USphcomp.box());
      cellQuotientValAndInvRadial(USphtemp,a_dVolr,dVolrInv,USphcomp,a_bx,a_blockBox);
      USphtemp.copyTo(JUcomp);
    }
  // cout << "In InterpBoundariesEuler: USph_i.box() = " << USph_i.box() << endl;
  // h5.writePatch(a_dxiPerp,a_JU,"USph:Block=" + to_string(a_block));
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void CubedSphereShell::
consToSphInterpEuler(                    
                        MBLevelBoxData<T,C,MEM,CTR>& a_JU,
                        MBInterpOp& a_iop,
                        MBLevelBoxData<T,1,MEM,CTR>& a_dVolr,
                        unsigned int a_order)
{
  //HDF5Handler h5;
  auto layout = a_JU.layout();
  Point size = layout.domain().getBlock(0).sizes();
  T dxiPerp = 1.0/size[1];
  for (auto dit : a_JU.layout())
    {
      unsigned int block = a_JU.layout().block(dit);
      auto& JU_i = a_JU[dit];
      auto& dVolr_i = a_dVolr[dit];
      Box bx_i = a_JU.layout()[dit];      
      //Box blockBox = a_JU.layout().getBlock(block).domain().box();
       Box blockBox = a_JU.layout()[dit];
      CubedSphereShell::
        consToSphNGEuler(JU_i,dVolr_i,bx_i,blockBox,dxiPerp,block,a_order);
      //h5.writePatch(dxiPerp, JU_i, "consToSphInterpEuler:USph:Block" + to_string(block));
    } 
  a_JU.exchange();
  // cout <<  "Building interpolation operator" << endl;
  // auto op = CubedSphereShell::InterpOp(a_JU,4);
  a_iop.apply(a_JU,a_JU);
}


template<typename T, unsigned int C, MemType MEM>
void CubedSphereShell::
consToCartNGEuler(                    
                 BoxData<T,C,MEM>& a_JU,
                 BoxData<T,1,MEM>& a_dVolr,
                 Box a_bx,
                 Box a_blockBox,
                 T a_dxiPerp,
                 int a_block,                 
                 unsigned int a_order)
{
  PR_TIMERS("consToCartNGEuler");
  
  //HDF5Handler h5;
  BoxData<T,C,MEM> UCart(a_bx);
  UCart.setVal(1.);
  // Get geometric data.
  Box bxMetric = a_bx.grow(1);
  BoxData<T,1,MEM> deta = detA<T,MEM>(bxMetric,a_dxiPerp,0); 
  BoxData<T,1,MEM> detaInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_deta)
     {
       a_output(0) = 1.0/a_deta(0);
     }
     ,deta);

  // Slice data and compute quotients.
  for (int comp = 0; comp < C; comp++)
  {
    BoxData<T,1,MEM> JCart = slice(a_JU,comp);
    BoxData<T,1,MEM> Cart = slice(UCart,comp);      
    cellQuotientValAndInv(Cart,deta,detaInv,JCart,a_bx,a_blockBox);
  }

  // scale out radial dVol.
  BoxData<T,1,MEM> dVolrInv = forall<T,1,MEM>
    ([ ] PROTO_LAMBDA
     (Var<T,1,MEM>& a_output, const Var<T,1,MEM>& a_dvolr)
     {
       a_output(0) = 1.0/a_dvolr(0);
     }
     ,a_dVolr); 
  for (int comp = 0; comp < C; comp++)
    {
      BoxData<T,1,MEM> UCartcomp = slice(UCart,comp);
      BoxData<T,1,MEM> JUcomp = slice(a_JU,comp);
      BoxData<T,1,MEM> UCarttemp(UCartcomp.box());
      cellQuotientValAndInvRadial(UCarttemp,a_dVolr,dVolrInv,UCartcomp,a_bx,a_blockBox);
      UCarttemp.copyTo(JUcomp);
    }
  // cout << "In InterpBoundariesEuler: USph_i.box() = " << USph_i.box() << endl;
  // h5.writePatch(a_dxiPerp,a_JU,"USph:Block=" + to_string(a_block));
}
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void CubedSphereShell::
consToCartInterpEuler(                    
                        MBLevelBoxData<T,C,MEM,CTR>& a_JU,
                        MBInterpOp& a_iop,
                        MBLevelBoxData<T,1,MEM,CTR>& a_dVolr,
                        unsigned int a_order)
{
  //HDF5Handler h5;
  auto layout = a_JU.layout();
  Point size = layout.domain().getBlock(0).sizes();
  T dxiPerp = 1.0/size[1];
  for (auto dit : a_JU.layout())
    {
      unsigned int block = a_JU.layout().block(dit);
      auto& JU_i = a_JU[dit];
      auto& dVolr_i = a_dVolr[dit];
      Box bx_i = a_JU.layout()[dit];      
      Box blockBox = a_JU.layout().getBlock(block).domain().box();
      CubedSphereShell::
        consToCartNGEuler(JU_i,dVolr_i,bx_i,blockBox,dxiPerp,block,a_order);
      //h5.writePatch(dxiPerp, JU_i, "consToCartInterpEuler:UCart:Block" + to_string(block));
    } 
  a_JU.exchange();
  // cout <<  "Building interpolation operator" << endl;
  // auto op = CubedSphereShell::InterpOp(a_JU,4);
  a_iop.apply(a_JU,a_JU);
}


#if 0
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void CubedSphereShell::
consToSphInterpEuler(                    
                        MBLevelBoxData<T,C,MEM,CTR>& a_JU,
                        MBLevelBoxData<T,1,MEM,CTR>& a_dVolr,
                        unsigned int a_order)
{
  auto iop = CubedSphereShell::InterpOp(a_JU,a_order);
  CubedSphereShell::consToSphInterpEuler(a_JU,iop,a_dVolr,a_dx,a_order);
}
#endif
template <typename T, MemType MEM>
void CubedSphereShell::
cellMatrixQuotientMatAndInv(
                            BoxData<T,DIM,MEM,1>& a_quot,
                            BoxData<T,DIM,MEM,1>& a_quotBar,
                            const BoxData<T,DIM,MEM,DIM>& a_denominator,
                            const BoxData<T,DIM,MEM,DIM>& a_denomInv,
                            const BoxData<T,DIM,MEM,1>& a_numerator)
{
  T one = 1.0;
  forallInPlace(f_matVecProd,a_quotBar,a_denomInv,a_numerator,one);
  BoxData<T,DIM,MEM,1> corr(a_quotBar.box().grow(-1));
  corr.setToZero();
  // Note that we do not compute correction in the zero direction.
  for (int dir = 1; dir < DIM; dir++)
    {
      BoxData<T,DIM,MEM,1> Dquotbar = Stencil<T>::Derivative(1,dir)(a_quotBar);
      BoxData<T,DIM,MEM,DIM> Ddenom = Stencil<T>::Derivative(1,dir)(a_denominator);
      BoxData<T,DIM,MEM,1> prod(Ddenom.box());
      T twelveth = -1.0/12.0;
      forallInPlace(f_matVecProd,prod,Ddenom,Dquotbar,twelveth);      
      corr += prod;
    }
  BoxData<T,DIM,MEM> temp2 = forall<T,DIM,MEM>
    ([ ] PROTO_LAMBDA(Var<T,DIM,MEM>& a_output,
                      const Var<T,DIM,MEM>& a_num,
                      const Var<T,DIM,MEM>& a_corr)
     {
       for (int dir = 0; dir < DIM; dir++)
         {
           a_output(dir) = a_num(dir) + a_corr(dir);
         }
     },a_numerator,corr);
  forallInPlace(f_matVecProd,a_quot,a_denomInv,temp2,one);
}
template <typename T, MemType MEM>
void CubedSphereShell::
cellMatrixQuotientMatAndInv(
                            BoxData<T,DIM,MEM,1>& a_quot,
                            const BoxData<T,DIM,MEM,DIM>& a_denominator,
                            const BoxData<T,DIM,MEM,DIM>& a_denomInv,
                            const BoxData<T,DIM,MEM,1>& a_numerator,
                            Box a_quotBox,
                            Box a_blockBox)
{
  T one = 1.0;
  BoxData<T,DIM,MEM,1> quotBar(a_quotBox);
  BoxData<T,DIM,MEM,1> corr(a_quotBox);
  forallInPlace(f_matVecProd,quotBar,a_denomInv,a_numerator,one);
  corr.setToZero();
  // Note that we do not compute correction in the zero direction.
  for (int dir = 1; dir < DIM; dir++)
    {
      BoxData<T,DIM,MEM,1> Dquotbar =
        CubedSphereShell::deriv3case(quotBar,a_quotBox,a_blockBox,dir);
      BoxData<T,DIM,MEM,DIM> Ddenom = Stencil<T>::Derivative(1,dir)(a_denominator);
      BoxData<T,DIM,MEM,1> prod(Ddenom.box());
      T twelveth = -1.0/12.0;
      forallInPlace(f_matVecProd,prod,Ddenom,Dquotbar,twelveth);      
      corr += prod;
    }
  BoxData<T,DIM,MEM> temp2 = forall<T,DIM,MEM>
    ([ ] PROTO_LAMBDA(Var<T,DIM,MEM>& a_output,
                      const Var<T,DIM,MEM>& a_num,
                      const Var<T,DIM,MEM>& a_corr)
     {
       for (int dir = 0; dir < DIM; dir++)
         {
           a_output(dir) = a_num(dir) + a_corr(dir);
         }
     },a_numerator,corr);
  forallInPlace(f_matVecProd,a_quot,a_denomInv,temp2,one);
}
template <typename T, MemType MEM>
void CubedSphereShell::cellQuotientValAndInv(
                                             BoxData<T,1,MEM>& a_quot,
                                             BoxData<T,1,MEM>& a_quotBar,

                                             const BoxData<T,1,MEM>& a_denominator,
                                             const BoxData<T,1,MEM>& a_denomInv,
                                             const BoxData<T,1,MEM>& a_numerator)
{
  a_quot.setToZero();
  a_quotBar.setToZero();
  auto temp = a_denomInv*a_numerator;
  temp.copyTo(a_quotBar);
  BoxData<T,1,MEM> corr(a_quotBar.box().grow(-1));
  corr.setToZero();
  T twelveth = -1.0/12.0;
  // Note that we do not compute correction in the zero direction.
  for (int dir = 1; dir < DIM; dir++)
    {
      BoxData<T,1,MEM> Dquotbar = Stencil<T>::Derivative(1,dir)(a_quotBar);
      BoxData<T,1,MEM> Ddenom = Stencil<T>::Derivative(1,dir)(a_denominator);
      BoxData<T,1,MEM> prod = Dquotbar*Ddenom;      
      corr += prod;
    }
  corr *= twelveth;
  a_quotBar.copyTo(a_quot);
  a_quot += corr;
}
template <typename T, MemType MEM>
void CubedSphereShell::cellQuotientValAndInv(
                                             BoxData<T,1,MEM>& a_quot,
                                             const BoxData<T,1,MEM>& a_denominator,
                                             const BoxData<T,1,MEM>& a_denomInv,
                                             const BoxData<T,1,MEM>& a_numerator,
                                             Box a_quotBox,
                                             Box a_blockBox)
{
  a_quot.setToZero();
  BoxData<T,1,MEM> quotBar = a_denomInv*a_numerator;
  BoxData<T,1,MEM> corr(a_quotBox);
  corr.setToZero();
  T twelveth = -1.0/12.0;
  for (int dir = 1; dir < DIM; dir++)
    {
      BoxData<T,1,MEM> Dquotbar =
        CubedSphereShell::deriv3case(quotBar,a_quotBox,a_blockBox,dir);
      BoxData<T,1,MEM> Ddenom = Stencil<T>::Derivative(1,dir)(a_denominator);
      BoxData<T,1,MEM> prod = Dquotbar*Ddenom;      
      corr += prod;
    }
  corr *= twelveth;
  quotBar.copyTo(a_quot);
  a_quot += corr;
}
template <typename T, MemType MEM>
void CubedSphereShell::
cellQuotientValAndInvRadial(
                            BoxData<T,1,MEM>& a_quot,
                            const BoxData<T,1,MEM>& a_denominator,
                            const BoxData<T,1,MEM>& a_denomInv,
                            const BoxData<T,1,MEM>& a_numerator,
                            Box a_quotBox,
                            Box a_blockBox)
{
  a_quot.setToZero();
  BoxData<T,1,MEM> quotBar = a_denomInv*a_numerator;
  T twelveth = -1.0/12.0;
  BoxData<T,1,MEM> Dquotbar =
    CubedSphereShell::deriv3case(quotBar,a_quotBox,a_blockBox,0);
  BoxData<T,1,MEM> Ddenom = Stencil<T>::Derivative(1,0)(a_denominator);
  BoxData<T,1,MEM> corr = Dquotbar*Ddenom;      
  corr *= twelveth;
  corr *= a_denomInv;
  quotBar.copyTo(a_quot);
  a_quot += corr;
}
template<typename T, MemType MEM>
inline BoxData<T,DIM,MEM,DIM>
CubedSphereShell::detA_A(
                         Box a_bx,
                         T a_h,
                         unsigned int a_block,
                         unsigned int a_centering) 
{
  PR_TIMERS("detA_A");
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights_0={1.,0.,0.,0.,0.};
  int qnum = 5;
  Box bxCentered = a_bx.extrude(Point::Basis(a_centering));
  int q1d1 = 5;
  int q1d2 = 5;
  Array<double,5> weight1 = weights;
  Array<double,5> weight2 = weights;
  if (a_centering == 1)
    {
      q1d1 = 2;
      weight1 = weights_0;
    }
  if (a_centering == 2)
    {
      q1d2 = 2;
      weight2 = weights_0;
    }
  Point low0 = a_bx.low();
  Point high0 = a_bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);  
  BoxData<T,DIM,MEM,DIM> output(a_bx);
  forallInPlace_p
    ([ ] PROTO_LAMBDA
     (Point a_pt,
      Var<T,DIM,MEM,DIM>& a_deta_amat,
      T a_h,
      Array<T,5> a_wgt1,
      Array<T,5> a_wgt2,
      int a_q1d1,
      int a_q1d2,
      Array<unsigned int,DIM> a_permute,
      Array<int,DIM> a_sign)
     {
       matType deta_amatav;
       for (int dir0 = 0; dir0 < DIM; dir0++)
         {         
           deta_amatav[dir0][0] = 0.;
           deta_amatav[dir0][1] = 0.;
           deta_amatav[dir0][2] = 0.;
         } 
       for (int k1 = 0; k1 < a_q1d1;k1++)
         {
           T offseta = (1.0*k1)/(a_q1d1-1);
           for (int k2 = 0; k2 < a_q1d2;k2++)
             {
               T offsetb = (1.0*k2)/(a_q1d2-1);
               matType deta_amat;
               
               f_detA_A_(a_pt,deta_amat,a_permute,a_sign,a_h,offseta,offsetb);
               for (int dir0=0; dir0 < DIM; dir0++)
                 {
                   for (int dir1=0; dir1 < DIM; dir1++)
                     {
                       deta_amatav[dir0][dir1] +=
                         deta_amat[dir0][dir1]*a_wgt1[k1]*a_wgt2[k2];
                     }
                 }
             }
         }
       for (int dir0 = 0; dir0 < DIM; dir0++)
         {
           a_deta_amat(dir0,0) = deta_amatav[dir0][0];
           a_deta_amat(dir0,1) = deta_amatav[dir0][1];
           a_deta_amat(dir0,2) = deta_amatav[dir0][2];          
         } 
       
     },bx0,output,a_h,weight1,weight2,q1d1,q1d2,permute[a_block],sign[a_block]);
  Operator::spreadSlice(output);
  if (a_centering == 0)
    {
      BoxData<T,DIM,MEM,1> outcol(a_bx);
      outcol.setToZero();
      Box bx2 = a_bx.extrude(Point::Basis(2));
      Box bx1 = a_bx.extrude(Point::Basis(1));
      BoxData<T,DIM,MEM> d1n_x_n(bx2);
      BoxData<T,DIM,MEM> n_x_d2n(bx1);
      nCrossD2n(n_x_d2n,a_bx,a_h,a_block);
      d1nCrossN(d1n_x_n,a_bx,a_h,a_block);
      T halfbyh = .5/a_h;
      outcol += Stencil<T>::FluxDivergence(1)(n_x_d2n,halfbyh);
      outcol += Stencil<T>::FluxDivergence(2)(d1n_x_n,halfbyh);
      forallInPlace([ ] PROTO_LAMBDA(
                                     Var<T,DIM,MEM,1>& a_d1n_x_d2n,
                                     Var<T,DIM,MEM,DIM>& a_deta_amat)
                    {
                      for (int dir = 0; dir < DIM; dir++)
                        {
                          a_deta_amat(dir,0) = a_d1n_x_d2n(dir);
                        }
                    },outcol,output);     
    }  
  return output;
}
#if 1
template<typename T, MemType MEM>
inline BoxData<T,DIM,MEM,DIM> CubedSphereShell::AMatrix(
                                                        Box a_bx,
                                                        T a_h,
                                                        unsigned int a_block,
                                                        unsigned int a_centering)
{
  PR_TIMERS("AMatrix");
  Array<Array<uint,DIM>,6>
    permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights_0={1.,0.,0.,0.,0.};
  int qnum = 5;
  Box bxCentered = a_bx.extrude(Point::Basis(a_centering));
  int q1d1 = 5;
  int q1d2 = 5;
  Array<double,5> weight1 = weights;
  Array<double,5> weight2 = weights;
  if (a_centering == 1)
    {
      q1d1 = 2;
      weight1 = weights_0;
    }
  if (a_centering == 2)
    {
      q1d2 = 2;
      weight2 = weights_0;
    }
  Point low0 = bxCentered.low();
  Point high0 = bxCentered.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  BoxData<T,DIM,MEM,DIM> output0 = forall_p<T,DIM,MEM,DIM>
    ([ ] PROTO_LAMBDA
     (Point a_pt,
      Var<T,DIM,MEM,DIM>& a_amat,
      T a_h,
      Array<T,5> a_wgt1,
      Array<T,5> a_wgt2,
      int a_q1d1,
      int a_q1d2,
      Array<unsigned int,DIM> a_permute,
      Array<int,DIM> a_sign)
     {
       for (int dir0 = 0; dir0 < DIM; dir0++)
         {
           a_amat(dir0,0) = 0.;
           a_amat(dir0,1) = 0.;
           a_amat(dir0,2) = 0.;
         } 
       for (int k1 = 0; k1 < a_q1d1;k1++)
         {
           T offseta = (1.0*k1)/(a_q1d1-1);
           for (int k2 = 0; k2 < a_q1d2;k2++)
             {
               T offsetb = (1.0*k2)/(a_q1d2-1);
               matType amat2;
               f_Amatrix_(a_pt,amat2,a_permute,a_sign,a_h,offseta,offsetb);
               for (int dir0=0; dir0 < DIM; dir0++)
                 {
                   for (int dir1=0; dir1 < DIM; dir1++)
                     {
                       a_amat(dir0,dir1) +=
                         amat2[dir0][dir1]*a_wgt1[k1]*a_wgt2[k2];
                     }
                 }
             }
         }
     },bx0,a_h,weight1,weight2,q1d1,q1d2,permute[a_block],sign[a_block]);
  BoxData<T,DIM,MEM,DIM> output(bxCentered);
  output0.copyTo(output);
  spreadSlice(output);
  return output;
}
#endif
template<typename T,MemType MEM>
BoxData<T,1,MEM> CubedSphereShell::detA(
                                        Box a_bx,
                                        T a_h,
                                        unsigned int a_centering)
{
  PR_TIMERS("detA");
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights_0={1.,0.,0.,0.,0.};
  int q1d1 = 5;
  int q1d2 = 5;
  Array<double,5> weight1 = weights;
  Array<double,5> weight2 = weights;
  if (a_centering == 1)
    {
      q1d1 = 2;
      weight1 = weights_0;
    }
  if (a_centering == 2)
    {
      q1d2 = 2;
      weight2 = weights_0;
     }
  Point low0 = a_bx.low();
  Point high0 = a_bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  BoxData<T,1,MEM> output0 =
    forall_p<T,1,MEM>([ ] PROTO_LAMBDA
                      (                  
                       Point a_pt,
                       Var<T,1,MEM>& a_out,
                       T a_h,
                       Array<T,5> a_wgt1,
                       Array<T,5> a_wgt2,
                       int a_q1d1,
                       int a_q1d2)
                      {
                        a_out(0) = 0.0;
                        for (int k1 = 0; k1 < a_q1d1; k1++)
                          {
                            T offseta = (1.0*k1)/(a_q1d1-1);
                            for (int k2 = 0; k2 < a_q1d2;k2++)
                              {
                                T offsetb = (1.0*k2)/(a_q1d2-1);
                                T outPoint;
                                f_detA_<T,MEM>(a_pt,outPoint,a_h,offseta,offsetb);  
                                a_out(0) += outPoint*a_wgt1[k1]*a_wgt2[k2]; 
                              }
                          }                        
                      },bx0,a_h,weight1,weight2,q1d1,q1d2);
  BoxData<T,1,MEM> output(a_bx);
  output0.copyTo(output);
  spreadSlice(output);
  return output;
}
template<typename T, MemType MEM>
BoxData<T,DIM,MEM> CubedSphereShell::adjA(
                                          Box a_bx,
                                          T a_h,
                                          int a_block,
                                          int a_centering)
{
  PR_TIMERS("adjA");
  BoxData<T,DIM,MEM,1> output(a_bx);
  T hInvHalf = .5/a_h;
  switch(a_centering)
    {
    case 0:
      {
        Box bx2 = a_bx.extrude(Point::Basis(2));
        Box bx1 = a_bx.extrude(Point::Basis(1));
        BoxData<T,DIM,MEM> d1n_x_n(bx2);
        BoxData<T,DIM,MEM> n_x_d2n(bx1);
        nCrossD2n(n_x_d2n,a_bx,a_h,a_block);
        d1nCrossN(d1n_x_n,a_bx,a_h,a_block);
        output |= Stencil<T>::FluxDivergence(1)(n_x_d2n,hInvHalf);
        output += Stencil<T>::FluxDivergence(2)(d1n_x_n,hInvHalf);
      }
      break;
    case 1:
      nCrossD2n(output,a_bx,a_h,a_block);
      output *= -1.0;
      break;
    case 2:
      d1nCrossN(output,a_bx,a_h,a_block);
      output *= -1.0;
      break;
    default:abort();
    }
  return output;
}
template<typename T, MemType MEM>
inline void
CubedSphereShell::d1nCrossN(
                            BoxData<T,DIM,MEM>& a_d1n_x_n,
                            Box a_bx,
                            T a_h,
                            uint a_block)
{
  PR_TIMERS("d1nCrossN");
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights0={1.0,0.,0.,0.,0.};
 Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  int qnum = 5;
  int qnum0 = 2;
  Box bx = a_d1n_x_n.box();
  Point low0 = bx.low();  
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  
  BoxData<T,DIM,MEM> tmp = 
    forall_p<T,DIM,MEM>
             ([]PROTO_LAMBDA
              (
               Point a_pt,
               Var<T,DIM,MEM>& a_d1n_x_n,
               Array<T,5>& a_weights,
               T a_h,
               Array<uint,DIM> a_permute,
               Array<int,DIM> a_sign,
               int a_qnum)
              {
                a_d1n_x_n(0) = 0.0;
                a_d1n_x_n(1) = 0.0;
                a_d1n_x_n(2) = 0.0;
                T beta = -.25*M_PI + .5*M_PI*a_pt[2]*a_h ;
                T tanb = tan(beta);
                Array<T,5> cross = {0.,0.,0.};
                for (int k = 0; k < a_qnum;k++)
                  {
                    T alpha = -.25*M_PI + .5*M_PI*(a_pt[1]*a_h + (k*a_h)/(a_qnum-1));
                    T tana = tan(alpha);
                    T secasq = 1.0 + tana*tana;
                    T normsq = M_PI*.5/(tana*tana + tanb*tanb + 1.0);       
                    cross[0] = secasq*tanb*normsq;
                    cross[1] = 0.;
                    cross[2] = -secasq*normsq;
                    a_d1n_x_n(0) += a_sign[0]*cross[a_permute[0]]*a_weights[k];
                    a_d1n_x_n(1) += a_sign[1]*cross[a_permute[1]]*a_weights[k];
                    a_d1n_x_n(2) += a_sign[2]*cross[a_permute[2]]*a_weights[k];
                  }       
              },bx0,weights,a_h,permute[a_block],sign[a_block],qnum);
   tmp.copyTo(a_d1n_x_n);
  // copy lowest slice to all the other slices. 
  Operator::spreadSlice(a_d1n_x_n);
}
template<typename T, MemType MEM>
inline void
CubedSphereShell::nCrossD2n(
                            BoxData<T,DIM,MEM>& a_n_x_d2n,
                            Box a_bx,
                            T a_h,
                            uint a_block)
{
  PR_TIMERS("nCrossD2n");
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights0={1.0,0.,0.,0.,0.};
   Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  int qnum = 5;
  int qnum0 = 2;
  Box bx = a_n_x_d2n.box();
  Point low0 = bx.low();  
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  
  BoxData<T,DIM,MEM> tmp =
    forall_p<T,DIM,MEM>([]PROTO_LAMBDA
             (
              Point a_pt,
              Var<T,DIM,MEM>& a_n_x_d2n,
              Array<T,5>& a_weights,
              T a_h,
              Array<uint,DIM> a_permute,
              Array<int,DIM> a_sign,
              int a_qnum)
             {
               a_n_x_d2n(0) = 0.0;
               a_n_x_d2n(1) = 0.0;
               a_n_x_d2n(2) = 0.0;
               T alpha = -.25*M_PI + .5*M_PI*a_pt[1]*a_h ;
               T tana = tan(alpha);
               Array<T,5> cross = {0.,0.,0.};
               for (int k = 0; k < a_qnum;k++)
                 {
                   T beta =
                     -.25*M_PI + .5*M_PI*(a_pt[2]*a_h+ (k*a_h)/(a_qnum-1));
                   T tanb = tan(beta);
                   T secbsq = 1.0 + tanb*tanb;
                   T normsq = M_PI*.5/(tana*tana + tanb*tanb + 1.0);       
                   cross[0] = secbsq*tana*normsq;
                   cross[1] = -secbsq*normsq;
                   cross[2] = 0.;
                   a_n_x_d2n(0) += a_sign[0]*cross[a_permute[0]]*a_weights[k];
                   a_n_x_d2n(1) += a_sign[1]*cross[a_permute[1]]*a_weights[k];
                   a_n_x_d2n(2) += a_sign[2]*cross[a_permute[2]]*a_weights[k];
                 }
             },bx0,weights,a_h,perm(a_block),permsign(a_block),qnum);
  tmp.copyTo(a_n_x_d2n);
  // copy lowest slice to all the other slices. 
  Operator::spreadSlice(a_n_x_d2n);
}
template <typename T, unsigned int C, MemType MEM>
BoxData<T,C,MEM>
CubedSphereShell::deriv3case(
                             const BoxData<T,C,MEM>& a_input,
                             Box a_outputBox,
                             Box a_validBox,
                             uint a_dir)
{
  PR_TIMERS("deriv3case");
  BoxData<T,C,MEM> output(a_outputBox);
  Stencil<T> loStencil = (-1.5)*Shift(Point::Zeros()) +
    2.0*Shift(Point::Basis(a_dir)) + (-.5)*Shift(Point::Basis(a_dir,2));
  Stencil<T> hiStencil = (1.5)*Shift(Point::Zeros()) +
    (-2.0)*Shift(Point::Basis(a_dir,-1)) + .5*Shift(Point::Basis(a_dir,-2));
  Box loGhostBox = a_outputBox.adjacent(a_dir,Side::Lo,1);
  Box hiGhostBox = a_outputBox.adjacent(a_dir,Side::Hi,1);  
  Box loBox = a_outputBox.face(a_dir,Side::Lo,1);
  Box hiBox = a_outputBox.face(a_dir,Side::Hi,1);
  output |= Stencil<T>::Derivative(1,a_dir,2)(a_input);
  if (!a_validBox.contains(loGhostBox))
    {     
      output |= loStencil(a_input,loBox);
    }
  if (!a_validBox.contains(hiGhostBox))
    {      
      output |= hiStencil(a_input,hiBox);
    }
    return output;
}
