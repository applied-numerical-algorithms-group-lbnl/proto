//==============================================================
// CUBED SPHERE SHELL MAP
template<MemType MEM>
void MBMap_CubedSphereShell<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT)
{
    auto& dx = this->dx();
    auto block = this->block();
#if 1
    auto C2C = Stencil<double>::CornersToCells(4);
    int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    double r0 = CUBED_SPHERE_SHELL_R0;
    double r1 = CUBED_SPHERE_SHELL_R1;
    Box b = C2C.domain(a_X.box());
    BoxData<double, DIM, MEM> X0(b);
    BoxData<double, 1, MEM> J0(b);
    forallInPlace_p(
            f_cubeSphereMap, X0, J0, dx, r0, r1, block, r_dir);
    X0.copyTo(a_X);
    a_J |= C2C(J0);
#else
    BoxData<double, 1, MEM> radius(a_X.box());
    forallInPlace_p(
            [] PROTO_LAMBDA (
                Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
            {
            a_R(0) = a_pt[0]*a_dx[0] + CUBED_SPHERE_SHELL_R0;
            }, radius, dx);
    Operator::cubedSphereGeometry(
            a_X, a_NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], block);
#endif
}

//==============================================================
// CUBED SPHERE SHELL POLAR MAP
template<MemType MEM>
void MBMap_CubedSphereShellPolar<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT)
{
    auto& dx = this->dx();
    auto block = this->block();
#if 1
    int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    double r0 = CUBED_SPHERE_SHELL_R0;
    double r1 = CUBED_SPHERE_SHELL_R1;
    BoxData<double, DIM, MEM> X0(a_X.box().grow(1));
    BoxData<double, 1, MEM> J0(X0.box());
    BoxData<double, DIM, MEM> X1(a_X.box().grow(1));
    BoxData<double, 1, MEM> J1(X1.box());
    forallInPlace_p(
            f_cubeSphereMap, X0, J0, dx, r0, r1, block, r_dir);
    int focalBlock = this->focalBlock();
    forallInPlace_p(
            f_cartesianToSphericalMap, X1, J1, X0, focalBlock, r_dir); 
    J0.copyTo(a_J);
    a_J *= J1;
    X1.copyTo(a_X);
#else

    BoxData<double, 1, MEM> radius(a_X.box());
    forallInPlace_p(
            [] PROTO_LAMBDA (
                Point& a_pt, Var<double, 1, MEM>& a_R, Array<double, DIM> a_dx)
            {
            a_R(0) = a_pt[0]*a_dx[0] + CUBED_SPHERE_SHELL_R0;
            }, radius, dx);
    Array<BoxData<double,DIM,MEM>, DIM> NT;
    Operator::cubedSphereGeometry(
            a_X, NT, a_J, radius, a_J.box().grow(Point::Ones()), dx[1], a_block);
#endif
}

//==============================================================
// CUBED SPHERE SHELL UTILITY
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
        const Array<Point, DIM+1>   a_ghost)
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

template<template<typename, MemType> class OpType,
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
        const Array<Point, DIM+1>& a_ghost,
        unsigned int a_order)
{
    std::vector<Point> footprint;
    Box K = Box::Kernel(2);
    for (auto pi : K)
    {
        // Assumes radial coordinate is the first
        if (pi.abs().sum() <= 2 
            && pi[CUBED_SPHERE_SHELL_RADIAL_COORD] == 0)
        {
            footprint.push_back(pi);
        }
    }
    MBInterpOp op(a_ghost[0], a_order);
    for (int bi = 0; bi < 6; bi++)
    {
        MBLevelMap<MBMap_CubedSphereShellPolar, MEM> polarMap;
        polarMap.define(a_layout, a_ghost, bi);
        op.define(polarMap, footprint, bi); 
    }
    return op; 
}
    
template<typename T, unsigned int C, MemType MEM, Centering CTR>
MBInterpOp
CubedSphereShell::InterpOp(
        MBLevelBoxData<T,C,MEM,CTR>& a_data,
        unsigned int a_order)
{
    return InterpOp<MEM>(a_data.layout(), a_data.ghost(), a_order);
}
    
template<typename T, unsigned int C, MemType MEM, Centering CTR>
void
CubedSphereShell::InterpBoundaries(
        MBLevelBoxData<T,C,MEM,CTR>& a_data,
        unsigned int a_order)
{
    auto op = InterpOp(a_data, a_order);
    op.apply(a_data, a_data);
}
template <typename T, unsigned int C, MemType MEM>
void cellMatrixQuotientMatAndInv(
                                 BoxData<T,C,MEM,1>& a_quot,
                                 BoxData<T,C,MEM,1>& a_quotBar,
                                 const BoxData<T,C,MEM,C>& a_denominator,
                                 const BoxData<T,C,MEM,C>& a_denomInv,
                                 const BoxData<T,C,MEM,1>& a_numerator)
{
  a_quotBar.define(a_denomInv.box());
  T one = 1.0;
  forallInPlace(f_matVecProd,a_quotBar,a_denomInv,a_numerator,one);
  BoxData<T,C,MEM,1> corr(a_quotBar.box().grow(-1));
  corr.setToZero();
  for (int dir = 0; dir < DIM; dir++)
    {
      BoxData<T,C,MEM,1> Dquotbar = Stencil<T>::Derivative(1,dir)(a_quotBar);
      BoxData<T,C,MEM,C> Ddenom = Stencil<T>::Derivative(1,dir)(a_denominator);
      BoxData<T,C,MEM,1> prod(Ddenom.box());
      T twelveth = -1.0/12.0;
      forallInPlace(f_matVecProd,prod,Ddenom,Dquotbar,twelveth);      
      corr += prod;
    }
  a_quot = a_quotBar + corr;
}
template<typename T, MemType MEM>
inline void detA_A(
                   BoxData<T,DIM,MEM,DIM>& a_detA_A,
                   Box a_bx,
                   T a_h,
                   unsigned int a_block,
                   unsigned int a_centering)
{
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights_0={1.,0.,0.,0.,0.};
  Array<Array<uint,DIM>,6>
     permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  int qnum = 5;
  Point low0 = a_bx.grow(2).low();
  low0[0] = a_bx.low()[0];
  Point high0 = a_bx.grow(2).high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  Box bx2 = a_bx.extrude(Point::Basis(2));
  Box bx1 = a_bx.extrude(Point::Basis(1));
  BoxData<T,DIM,MEM> d1n_x_n(bx2);
  BoxData<T,DIM,MEM> n_x_d2n(bx1);
  Stencil<T> stl1 = Stencil<T>::FluxDivergence(1);
  Stencil<T> stl2 = Stencil<T>::FluxDivergence(2); 
  Box bxdef = stl1.range(n_x_d2n.box()) & stl2.range(d1n_x_n.box());
  BoxData<T,DIM,MEM> d1n_x_d2n(a_bx);
  T dxinv = .5/a_h;
  d1n_x_d2n += stl1(n_x_d2n,dxinv);
  d1n_x_d2n += stl2(d1n_x_n,dxinv);
  BoxData<T,DIM,MEM,DIM> detA_ASlice(a_bx);
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
  forallInPlace_p([ ] PROTO_LAMBDA
                  (Point a_pt,
                   Var<T,DIM,MEM,DIM>& a_deta_aslice,
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
                        a_deta_aslice(dir0,0) = 0.;
                        a_deta_aslice(dir0,1) = 0.;
                        a_deta_aslice(dir0,2) = 0.;
                      } 
                    for (int k1 = 0; k1 < a_q1d1;k1++)
                      {
                        T offseta = (1.0*k1)/(a_q1d1-1);
                        for (int k2 = 0; k2 < a_q1d2;k2++)
                          {
                            T offsetb = (1.0*k2)/(a_q1d2-1);
                            //Var<T,DIM,MEM,DIM> deta_amat; 
                            matType deta_amat;
                            f_detA_A_(a_pt,deta_amat,a_h,offseta,offsetb);
                            for (int dir0=0; dir0 < DIM; dir0++)
                              {
                                for (int dir1=0; dir1 < DIM; dir1++)
                                  {
                                    a_deta_aslice(dir0,dir1) +=
                                      a_sign[dir0]*deta_amat[a_permute[dir0]][dir1]
                                      *a_wgt1[k1]*a_wgt2[k2];
                                  }
                              }
                          }
                      }
                  },a_bx,detA_ASlice,a_h,weight1,weight2,q1d1,q1d2,
                  permute[a_block],sign[a_block]);
  d1nCrossN(d1n_x_n,a_bx,a_h,a_block);
  nCrossD2n(n_x_d2n,a_bx,a_h,a_block);
  d1n_x_d2n.setToZero();
  d1n_x_d2n += stl1(n_x_d2n,dxinv);
  d1n_x_d2n += stl2(d1n_x_n,dxinv);
  
  int lowrow = 0;
  if (a_centering == 0)
    {
      lowrow = 1;
      for (int dir0 = 0; dir0 < DIM;dir0++)
        {
          auto slcDest = slice(a_detA_A,dir0,0,0);
          slcDest.setToZero(); 
          auto slcSrc = slice(d1n_x_d2n,dir0,0,0);          
          slcSrc.copyTo(slcDest);
        }
    }
  //  Operator::spreadSlice(a_detA_A);
  for (int dir0 = 0; dir0 < DIM;dir0++)
    {
      for (int dir1 = lowrow; dir1 < DIM; dir1++)
        {
          auto slcDest = slice(a_detA_A,dir0,dir1,0);
          auto slcSrc = slice(detA_ASlice,dir0,dir1,0);
          slcSrc.copyTo(slcDest);          
        }
    }
  // Operator::spreadSlice(a_detA_A);
}
template<typename T,MemType MEM>
inline void
scaledDetA_A(
              BoxData<T,DIM,MEM,DIM>& a_scaledDetA_A,
              const BoxData<T,1,MEM>& a_radialScale,
              Box a_bx,
              T a_h,
              unsigned int a_block,
              unsigned int a_centering,
              bool a_scaled)
{
 Array<Array<uint,DIM>,6>
    permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{-1,1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  detA_A(a_scaledDetA_A,a_bx,a_h,a_block,a_centering);
  return;
  if (a_scaled)
    {
      forallInPlace([] PROTO_LAMBDA
                    (Var<T,DIM,MEM,DIM>& a_scaleddeta_a,
                     const Var<T>& a_radialscale)
                    {
                      for (int dir0 = 0; dir0 < DIM; dir0++)
                        {
                          for (int dir1 = 0; dir1 < DIM; dir1++)
                            {
                              a_scaleddeta_a(dir0,dir1) *= a_radialscale(0);
                            }
                        }
                    },a_scaledDetA_A,a_radialScale);
    }
}
template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
inline void 
scaledDetA(
           BoxData<T,C,MEM,D,E>& a_output,
           const BoxData<T,1,MEM>& a_scale,
           Box a_bx,
           T a_h,
           unsigned int a_block,
           unsigned int a_centering,
           bool a_useScaling)
{
  Array<Array<uint,DIM>,6>
    permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{-1,1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  Box bxCentered;
  if (a_centering == 0)
    {
      bxCentered = a_bx.extrude(Point::Basis(0));
    } else
    {
      bxCentered = a_bx.extrude(Point::Basis(0) + Point::Basis(a_centering));
    }
  Point low0 = bxCentered.low();
  low0[0] = a_bx.low()[0];
  Point high0 = bxCentered.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  // a_output.define(bxCentered);
  BoxData<T,C,MEM,D,E> outTemp(bxCentered);  
  quadrature_detA(outTemp,bx0,a_h,a_centering);
  Operator::spreadSlice(outTemp);
  //forallInPlace_p([ ] PROTO_LAMBDA(Point a_pt,Var<T,1,MEM>& a_out){a_out(0) = a_pt[1];},bxCentered,outTemp);
  outTemp.copyTo(a_output);
  if (!a_useScaling) return;
  for (int ee = 0; ee < E; ee++)
      {
        for (int dd = 0; dd < D; dd++)
          {
            for (int cc = 0; cc < C; cc++)
              {
                BoxData<T> outTempSlice = slice(outTemp,cc,dd,ee);
                if (a_useScaling)
                  {
                    outTempSlice *= a_scale;
                  }
                if (C == DIM)
                  {
                    outTempSlice *= (T (sign[a_block][cc]));
                    auto outputSlice = slice(a_output,permute[a_block][cc],dd,ee);
                    outTempSlice.copyTo(outputSlice);
                  }                 
              }
          }
      }
}
template<typename T,MemType MEM>
inline void 
quadrature_detA(
                BoxData<T,1,MEM>& a_output,
                Box a_bx,
                T a_h,
                unsigned int a_centering)
{
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights_0={1.,0.,0.,0.,0.};
  int qnum = 5;  
  int q1d1 = qnum;
  int q1d2 = qnum;
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
  forallInPlace_p([ ] PROTO_LAMBDA
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
                            f_detA_(a_pt,outPoint,a_h,offseta,offsetb);  
                            a_out(0) += outPoint*a_wgt1[k1]*a_wgt2[k2]; 
                          }
                    }
                  },a_bx,a_output,a_h,weight1,weight2,q1d1,q1d2);
}
template<typename T,MemType MEM>
inline void NTRowCS(
                    BoxData<T,DIM,MEM>& a_NTRow,
                    const BoxData<T,1,MEM>& a_adjD,
                    Box a_bx,
                    T a_h,
                    int a_block,
                    int a_centering)
{
  a_NTRow.define(a_bx);
  T half = .5;
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
        a_NTRow |= Stencil<T>::FluxDivergence(1)(n_x_d2n,half);
        a_NTRow += Stencil<T>::FluxDivergence(2)(d1n_x_n,half);
      }
      break;
    case 1:
      nCrossD2n(a_NTRow,a_bx,a_h,a_block);
      break;
    case 2:
      d1nCrossN(a_NTRow,a_bx,a_h,a_block);
      break;
    default:abort();
    }
  forallInPlace([] PROTO_LAMBDA(
                                Var<T,DIM>& a_ntrow,
                                const Var<T>& a_adjd)
                {
                  for (int dir = 0; dir < DIM; dir++)
                    {
                      a_ntrow(dir) *= a_adjd(0);
                    }
                },a_NTRow,a_adjD);
}
template <typename T, MemType MEM>
inline void d1nCrossN(
               BoxData<T,DIM,MEM>& a_d1n_x_n,
               Box a_bx,
               T a_h,
               uint a_block)
{
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
  
  forallInPlace_p
    ([]PROTO_LAMBDA
     (Point a_pt,
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
           T alpha = -.25*M_PI + .5*M_PI*(a_pt[1]*a_h + k*a_h/(a_qnum-1));
           T tana = tan(alpha);
           T secasq =1.0 + tan(alpha)*tan(alpha);
           T normsq = 1.0/(tana*tana + tanb*tanb + 1.0);       
           cross[0] = secasq*tanb*normsq;
           cross[1] = 0.;
           cross[2] = -secasq*normsq;
           a_d1n_x_n(0) += a_sign[0]*cross[a_permute[0]]*a_weights[k];
           a_d1n_x_n(1) += a_sign[1]*cross[a_permute[1]]*a_weights[k];
           a_d1n_x_n(2) += a_sign[2]*cross[a_permute[2]]*a_weights[k];          
         }       
     },bx,a_d1n_x_n,weights,a_h,permute[a_block],sign[a_block],qnum);
  // copy lowest slice to all the other slices. 
  // Operator::spreadSlice(a_d1n_x_n);
}
template<typename T, MemType MEM>
inline void nCrossD2n(
                      BoxData<T,DIM,MEM>& a_n_x_d2n,
                      Box a_bx,
                      T a_h,
                      uint a_block)
{
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
  
  forallInPlace_p
    ([]PROTO_LAMBDA
     (Point a_pt,
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
             -.25*M_PI + .5*M_PI*(a_pt[2]*a_h + k*a_h/(a_qnum-1));
           T tanb = tan(beta);
           T secbsq = 1.0 + tan(beta)*tan(beta);
           T normsq = 1.0/(tana*tana + tanb*tanb + 1.0);       
           cross[0] = secbsq*tana*normsq;
           cross[1] = -secbsq*normsq;
           cross[2] = 0.;
           a_n_x_d2n(0) += a_sign[0]*cross[a_permute[0]]*a_weights[k];
           a_n_x_d2n(1) += a_sign[1]*cross[a_permute[1]]*a_weights[k];
           a_n_x_d2n(2) += a_sign[2]*cross[a_permute[2]]*a_weights[k];
         }
     },bx,a_n_x_d2n,weights,a_h,permute[a_block],sign[a_block],qnum);
  // copy lowest slice to all the other slices. 
  // Operator::spreadSlice(a_n_x_d2n);
}
#if 0
template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
inline void 
quadratureTemplate(
                   BoxData<T,C,MEM,D,E>& a_output,
                   Box a_bx,
                   T a_h,
                   unsigned int a_centering)
{
  Array<T,5>
    weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
  Array<T,5>
    weights_0={1.,0.,0.,0.,0.};
  int qnum = 5;  
  int q1d1 = qnum;
  int q1d2 = qnum;
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
  forallInPlace_p([] PROTO_LAMBDA
                  (                  
                   Point a_pt,
                   Var<T,C,MEM,D,E>& a_out,
                   T a_h,
                   Array<T,5> a_wgt1,
                   Array<T,5> a_wgt2,
                   int a_q1d1,
                   int a_q1d2)
                  {                       
                    for (int k1 = 0; k1 < a_q1d1;k1++)
                      {
                        T offseta = (1.0*k1)/(a_q1d1-1);
                        for (int k2 = 0; k2 < a_q1d2;k2++)
                          {
                            T offsetb = (1.0*k2)/(a_q1d2-1);
                            Var<T,C,MEM,D,E> outPoint;
                            FUNC(a_pt,outPoint,a_h,offseta,offsetb);
                            for (int ee = 0 ; ee < E; ee++)
                              for (int dd = 0 ; dd < D; dd++)
                                for (int cc = 0 ; cc < C; cc++)
                                  a_out(cc,dd,ee) +=
                                    outPoint(cc,dd,ee)*a_wgt1[k1]*a_wgt2[k2];          
                          }
                      }
                  },a_bx,a_output,a_h,weight1,weight2,q1d1,q1d2);
}

template<typename T,MemType MEM,unsigned int C, unsigned int D,unsigned int E>
inline void 
scaled2DTo3D_template(
                      BoxData<T,C,MEM,D,E>& a_output,
                      BoxData<T,1,MEM,1,1>& a_scale,
                      Box a_bx,
                      T a_h,
                      unsigned int a_block,
                      unsigned int a_centering,
                      bool a_useScaling);
{
  Array<Array<uint,DIM>,6>
    permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
  Array<Array<int,DIM>,6>
    sign = {{-1,1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
  bxCentered = a_bx.extrude(Point::Basis(a_centering));
  Point low0 = bxCentered.low();
  low0[0] = a_bx.low()[0];
  Point high0 =bxCentered.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);    
  a_output.define(bxCentered);
  BoxData<T,C,MEM,D,E> outTemp(bxCentered);
  quadratureTemplate(outTemp,bx0,a_h,a_centering);
  Operator::spreadSlice(outTemp);   
  for (int ee = 0; ee < E; ee++)
    {
      for (int dd = 0; dd < D; dd++)
        {
          for (int cc = 0; cc < C; ee++)
            {
              outTempSlice = slice(outTemp,cc,dd,ee);
              if (a_useScaling)
                {
                  outTempSlice *= a_scale;
                }
              if (C == DIM)
                {
                  outTempSlice *= (T sign[a_block][cc]);
                  outTempSlice.copyTo(slice(a_output,permute[a_block][cc],dd,ee));
                }                 
            }
        }
    } 
}
#endif
