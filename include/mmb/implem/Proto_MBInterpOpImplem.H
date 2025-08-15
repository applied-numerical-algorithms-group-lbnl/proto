
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_binomialPower_tmp(Var<T, 1, MEM>& a_xp, const Var<T, DIM, MEM>& a_x, Point a_p)
{
    a_xp(0) = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (a_p[dir] != 0)
        {
            a_xp(0) *= pow(a_x(dir), a_p[dir]);
        }
    }
}
PROTO_KERNEL_END(f_binomialPower_tmp, f_binomialPower)

MBPointInterpOp::MBPointInterpOp(
        MBDataPoint              a_dst,
        std::vector<MBDataPoint> a_srcs,
        std::vector<double>      a_coefs)
{
    PROTO_ASSERT(a_srcs.size() == a_coefs.size(),
            "MBPointInterpOp::constructor | Error: source and coefs should be the same size.");
    m_dst = a_dst;
    m_srcs = a_srcs;
    m_coefs = a_coefs;
}

template<template<MemType> typename MAP, MemType MEM>
MBPointInterpOp::MBPointInterpOp(
                                 MBDataPoint a_dst,
                                 Point a_ghost,
                                 const MBLevelMap<MAP, MEM>& a_map,
                                 unsigned int a_order)
{
    PR_TIME("MBPointInterpOp::constructor");
    m_order = a_order;
    m_dst = a_dst;
    // cout << "constructor: m_dst  = " << m_dst.point() << " , " << m_dst.block() << endl;
    computeSourcePoints(a_map,a_ghost);
    auto layout = a_map.layout();
    computeExponents();    
    //validateMatrixDimensions(m_footprint);
    auto momentData = computeMoments(a_map);
    computeCoefficients(momentData,layout);
}

template<template<MemType> typename MAP, MemType MEM>
void MBPointInterpOp::computeSourcePoints(const MBLevelMap<MAP, MEM>& a_map,
                                          Point a_ghost)
{
  PR_TIME("MBPointInterpOp::computeSourcePoints");
  // compute footprint.  
  auto layout = a_map.layout();  
  auto domain = layout.domain();
  auto destIndex = m_dst.index();
  Box bx=layout[destIndex].grow(a_ghost);
  unsigned int destBlock = layout.block(destIndex);
  Box destBlockBox = layout.domain().getBlock(destBlock).box();
  Box destPatchBox = layout[destIndex].grow(a_ghost);
  Point lowvalid = destPatchBox.low();
  Point highvalid = destPatchBox.high();
  lowvalid[0] = max(lowvalid[0],0);
  highvalid[0] = min(highvalid[0],destBlockBox.high()[0]);
  destPatchBox = Box(lowvalid,highvalid);
  Box destPatchBoxSmall = layout[destIndex].grow(a_ghost - Point::Ones(2));
  lowvalid = destPatchBoxSmall.low();
  highvalid = destPatchBoxSmall.high();
  lowvalid[0] = max(lowvalid[0],0);
  highvalid[0] = min(highvalid[0],destBlockBox.high()[0]);
  destPatchBoxSmall = Box(lowvalid,highvalid);
  // cout << "m_dst  = " << m_dst.point() << " , " << m_dst.block() << endl;
  
  
  // Special for cubed-sphere, where the stencil is in the spherical direction only
  vector<Point> footprint0;
  footprint0.clear();
  Box kernel2d = Box(Point(0,-1,-1),Point(0,1,1));
  if (m_order == 4)
    {
      for (auto bit : kernel2d)
        {
          footprint0.push_back(bit);
        }
      footprint0.push_back(Point::Basis(1,2));
      footprint0.push_back(Point::Basis(1,-2));
      footprint0.push_back(Point::Basis(2,2));
      footprint0.push_back(Point::Basis(2,-2));
    }
  else
    {
      if (procID()==0) cout << "m_order != 4 - abort" << endl;
      abort();
    }
  MBPoint mbptTmp(m_dst.point(),m_dst.block());
  auto mbptCenter = a_map.op(0).findValid(domain,mbptTmp); 
  auto centerPoint = mbptCenter.point;
  auto centerBlock = mbptCenter.block;
  auto centerIndex = layout.find(centerPoint/layout.boxSize(centerBlock),centerBlock);
  auto centerBlockBox = layout.domain().getBlock(centerBlock).box();
  std::set<Point> footprint;
  footprint.clear();
  m_footprint.clear();
  m_srcs.clear();
  int count = -1;
  
  Point centerPullbackPoint = domain.
     convertPoint(centerPoint,centerBlock,m_dst.block(),PR_CELL);
  // if ((m_dst.point() == Point(0,5,36)) && (m_dst.block() == 2))
  //       {
  //         cout << "m_dst = " << m_dst << endl;
  //         cout << "centerPoint = " << centerPoint
  //              << ", centerBlock = " << centerBlock
  //              << ", centerPullbackPoint = " << centerPullbackPoint
  //              << ", centerIndex = " << centerIndex
  //              << endl;
  //       }
  for (auto fit : footprint0)
    {
      int srcBlock; 
      Point ptPullback;
      Point pt = fit + centerPoint;
      Point ptsrc;
      MBIndex srcIndex;
      Point boundaryDir;
      Point dir(0,0,0);
      if (centerBlockBox.containsPoint(pt))       
        {
          ptsrc = pt;
          ptPullback = domain.
            convertPoint(pt,centerBlock,m_dst.block(),PR_CELL);
          srcBlock = centerBlock;
          srcIndex = centerIndex;
        }
      else
        {
          //dir = layout[centerIndex].whichBoundaryContains(pt);
          dir = centerBlockBox.whichBoundaryContains(pt);
          if (dir.codim() == 1)
            {
              srcIndex =
              layout.findAdjacentIndex(layout.point(centerIndex),centerBlock,dir);
              srcBlock =layout.block(srcIndex);
              ptsrc = domain.convertPoint(pt,centerBlock,srcBlock,PR_CELL);
              ptPullback = domain.
                convertPoint(ptsrc,srcBlock,m_dst.block(),PR_CELL);
            }
          else
            {
              MBPoint mbpt(pt,centerBlock);
              auto mbptSrc = a_map.op(0).findValid(domain,mbpt);
              ptsrc = mbptSrc.point;
              srcBlock = mbptSrc.block;
              srcIndex = layout.
                find(ptsrc/layout.boxSize(srcBlock),srcBlock);
              ptPullback = domain.
                convertPoint(ptsrc,srcBlock,m_dst.block(),PR_CELL);
            }
        }
      if (!footprint.contains(ptPullback))
        {
          if (srcBlock == destBlock)
            {
              boundaryDir = Point::Zeros();
              MBDataPoint mbdpt(destIndex,ptPullback,layout,boundaryDir,destBlock);
              m_srcs.push_back(mbdpt);
              m_footprint.push_back(ptsrc);
              footprint.insert(ptPullback);
            }
          else
            {
              boundaryDir = layout[destIndex].whichBoundaryContains(ptPullback);
              MBDataPoint mbdpt(destIndex,ptPullback,layout,boundaryDir,srcBlock);
              m_srcs.push_back(mbdpt);
              m_footprint.push_back(ptsrc);
              footprint.insert(ptPullback);
            }
          // if ((m_dst.point() == Point(0,5,36)) && (m_dst.block() == 2))
          //   {
          //     count++;
          //     cout << "================= \n"
          //          << "count = " << count << "\n"
          //          << "srcIndex = " <<  srcIndex << "\n"
          //          << "destIndex = " << destIndex << "\n"
          //          << "ptsrc = " << ptsrc << "\n"
          //          << "boundaryDir = " << boundaryDir << "\n"
          //          << "pt = " << pt << "\n"
          //          << "ptPullback = " << ptPullback << "\n"
          //          << "src, dest blocks = " << srcBlock << ", " << destBlock << "\n"
          //          <<  "================= " << endl << endl;
          //   }  
        }
    }
    
      // else
      //   {
      //     cout << "duplicate point " << ptPullback << endl;
      //   }       
  // if (m_footprint.size() != footprint0.size())
  //   {
  //     cout << "center point = " << centerPoint
  //          << ", center block = " << centerBlock
  //          << ", destination block = " << m_dst.block()
  //          << ", footprint size = " << m_footprint.size()
  //          << endl;
  //   }
  //cout << endl;
}
void MBPointInterpOp::computeExponents()
{
    PR_TIME("MBPointInterpOp::computeExponents");
    Box sourceSpanBox;
    for (auto fi : m_footprint)
    {
      sourceSpanBox &= fi;
    }
    Point maxExponents = Point::Ones(m_order);

    // If the span of the footprint only contains a single layer
    // of data along the ith axis, assume that the interpolating
    // polynomial is in a plane orthogonal to that axis
    for (int ii = 0; ii < DIM; ii++)
    {
        if (sourceSpanBox.size(ii) == 1) { maxExponents[ii] = 1; }
    }
    m_exponents.clear();
    for (auto bi : Box(maxExponents))
    {
        if (bi.sum() < m_order)
        {
            m_exponents.push_back(bi);
        }
    }
}
void MBPointInterpOp::
validateMatrixDimensions(const std::vector<Point>& footprint) const
{
    PR_TIME("MBPointInterpOp::validateMatrixDimensions");
    int N = m_exponents.size();
    int M = m_srcs.size();
    if (N == 0)
    {
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Invalid system detected" << std::endl;
        pr_out() << "\tNo exponents were computed" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << m_dst << std::endl;
        pr_out() << "====================================================================" << std::endl;
    }
    if (M == 0)
    {
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Invalid system detected" << std::endl;
        pr_out() << "\tNo source points were found" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << m_dst << std::endl;
        pr_out() << "====================================================================" << std::endl;
    }
    if (N > M)
    {
        pr_out() << "====================================================================" << std::endl;
        pr_out() << "Underdetermined system detected" << std::endl;
        pr_out() << "\tN = " << N << ", M = " << M << std::endl;
        pr_out() << "\tDst Point: " << m_dst << std::endl;
        //pr_out() << "\tGhost: " << a_ghost << std::endl;
#if PR_VERBOSE > 1
        pr_out() << "\tFootprint: " << std::endl;
        for (auto fi : footprint)
        {
            pr_out() << "\t\t" << fi << std::endl;
        }
        pr_out() << "\tSrc Points: " << std::endl;
        for (auto si : m_srcs)
        {
            pr_out() << "\t\t" << si << std::endl;
        }
        pr_out() << "\tExponents: " << std::endl;
        for (auto ei : m_exponents)
        {
            pr_out() << "\t\t" << ei << std::endl;
        }
#endif
        pr_out() << "====================================================================" << std::endl;
    }
}

template<template<MemType> typename MAP, MemType MEM>
std::unordered_map<BlockIndex, std::vector<BoxData<double>>>
MBPointInterpOp::computeMoments(
    const MBLevelMap<MAP, MEM>& map) const
{
    PR_TIME("MBPointInterpOp::computeMoments");
    // the radius in physical space of the neighborhood containing
    // the source data points
    auto xCenter = map.apply(m_dst);
    double sourceDataRadius = computeSourceDataRadius(map, xCenter);

    std::unordered_map<BlockIndex, Box> momentDomains;
    std::unordered_map<BlockIndex, std::vector<BoxData<double>>> momentData;
    
    for (uint ii = 0; ii < m_footprint.size(); ii++)
    {
        int bi = m_srcs[ii].srcBlock();
        auto si = m_srcs[ii];
        if (momentDomains.find(bi) == momentDomains.end())
        {
          //momentDomains[bi] = Box(si.point(), si.point());
          momentDomains[bi] = Box(m_footprint[ii],m_footprint[ii]);
        } else {
            momentDomains[bi] &= m_footprint[ii];
        }
    }
    uint blockdest = m_dst.block();
    if (momentDomains.find(blockdest) == momentDomains.end())
      {
        momentDomains[blockdest] = Box(m_dst.point(),m_dst.point());
      }
    else
      {
        momentDomains[blockdest] &= m_dst.point();
      }
    auto C2C = Stencil<double>::CornersToCells(4);
    for (auto& [block, momentDomain] : momentDomains)
    {
        BoxData<double, DIM, MEM> Xi(C2C.domain(momentDomain).grow(1));
        // cout << "block = " << block << ", Box = "  << momentDomain << endl;
        BoxData<double, 1, MEM> Ji(Xi.box()); //Ji is unused
        map.apply(Xi, Ji, block);
        Xi -= xCenter;
        Xi /= sourceDataRadius;
        for (auto ei : m_exponents)
        {
            BoxData<double> moment = forall<double, 1>(f_binomialPower, Xi, ei);
            BoxData<double> momentAvg = C2C(moment);
            momentData[block].push_back(std::move(momentAvg));
        }
    }  
    return momentData;
}

template<template<MemType> typename MAP, MemType MEM>
double MBPointInterpOp::computeSourceDataRadius(
    const MBLevelMap<MAP, MEM>& map,
    Array<double, DIM>& xCenter) const
{
    PR_TIME("MBPointInterpOp::computeSourceDataRadius");
    double radius = 0.0;
    int counter = -1;
    for (auto& src : m_srcs)
    {
      counter++;
      auto xSrc = map.apply(src);
      auto dist = xSrc - xCenter;
      radius += dist.norm();
      if (isinf(radius))
        {
          cout << "radius = " << radius <<
            ", m_srcs.size() = " << m_srcs.size() <<
            ", dist = " << dist <<
            endl;
          cout << 
            ", counter = " << counter <<
            ", source point, block = " << src.point() << ", " <<
            src.block() <<
            ", dest point, block = " << m_dst.point() << ", " <<
            m_dst.block() << 
            endl;
          abort();
        }
    }
    radius /= (m_srcs.size());
    return radius;
}

void MBPointInterpOp::
computeCoefficients(
                    std::unordered_map
                    <BlockIndex, std::vector<BoxData<double>>>& momentData,
                    MBDisjointBoxLayout a_layout)
{
  PR_TIME("MBPointInterpOp::computeCoefficients");
  int N = m_exponents.size();
  int M = m_srcs.size();
  //cout << "N,M = " <<  N << ", " << M << endl;
  auto domain = a_layout.domain();
  Matrix<double> C(M,N);
  auto index = m_dst.index();
  for (int ii = 0; ii < M; ii++)
    {
      PR_TIMERS("build_C_matrix");
      auto& src = m_srcs[ii];
      auto srcBlock = src.srcBlock();
      for (int jj = 0; jj < N; jj++)
        {
          // C(ii,jj) = momentData[srcBlock][jj](src.point());
          C(ii,jj) = momentData[srcBlock][jj](m_footprint[ii]);
        }
    }
  Matrix<double> D(1,N);
  for (int jj = 0; jj < N; jj++)
    {
      auto destBlock = m_dst.block();
      PR_TIMERS("build_D_matrix");
      D(0,jj) = momentData[destBlock][jj](m_dst.point());
    }
  auto Cinv = C.inverse();
  auto S = D*Cinv;
  m_coefs.resize(M);
  double sum = 0.;
  for (int ii = 0; ii < M; ii++)
    {
      
      m_coefs[ii] = S(0,ii);
      sum += m_coefs[ii];
      // MBDataPoint src = m_srcs[ii];
      //     Point ptsrc = src.point();
      //     Point ptPullback =
      //       domain.convertPoint(ptsrc,src.block(),m_dst.block(),PR_CELL);
      //     Box bxDest = a_layout[m_dst.index()];
      //     auto dir = bxDest.whichBoundaryContains(ptPullback);
      //     auto mbdptsrc = MBDataPoint(m_dst.index(),ptPullback,a_layout,dir,src.block());
      //     m_srcs[ii] = mbdptsrc;
        
      //      //cout << "ii = "<< ii << ", m_srcs[ii] = " << m_srcs[ii] << endl;
    }
  if (abs(sum - 1.0) > 1.e-12)
  {
    pr_out() << "max(|sum - 1|) = " << abs(sum - 1.0) << " at " << m_dst.point() << ", " << m_dst.block() << ", size = " << M <<  endl;
  }
}

template<typename T, unsigned int C, MemType MEM>
void MBPointInterpOp::apply(
        MBLevelBoxData<T, C, MEM>& a_dst,
        MBLevelBoxData<T, C, MEM>& a_src)
{
  // PR_TIME("MBPointInterpOp::apply");
  // pr_out() << m_dst << endl;
  
  //auto dst = a_dst[m_dst]; //dst is a Var<T, C, MEM>
  MBIndex index = m_dst.index();
  Point pt = m_dst.point();
  auto dst = a_dst[index].var(pt);
  
  if ((m_dst.point() == Point(0,5,36)) && (m_dst.block() == 2))
      {
        // cout << "dest point = " << m_dst.point() << endl;
        auto layout = a_dst.layout();
        
        // T dstdummy = 0.;
        // for (uint ii = 0; ii < size(); ii++)
        //   {
        //     cout << "ii = " << ii
        //          << ", m_srcs[ii] = " << m_srcs[ii] << endl;
        //          cout << "weight = " << m_coefs[ii]
        //          << ", val = " << a_src[m_srcs[ii]](0)
        //          << endl;
        //     dstdummy += m_coefs[ii]*a_src[m_srcs[ii]](0);
        //   }
        // cout << "output value = " << dstdummy << endl;
        auto& bounddata = a_src.bound(m_srcs[0].index(),m_srcs[0].boundaryDir(),
             m_srcs[0].srcBlock());
        Box bx = layout[index];
        uint size = bx.high()[2]+ 1 - bx.low()[2];
        uint global = index.global();
        HDF5Handler h5;
        h5.writePatch(1.0,*bounddata.localData,"srcBoundData_"
                      + to_string(size) + "_" + to_string(global));
      }
    for (int ii = 0; ii < size(); ii++)
    {
        auto src = a_src[m_srcs[ii]]; //src is a Var<T, C, MEM>
        for (int cc = 0; cc < C; cc++)
        {
            if (ii == 0) { dst(cc) = 0.0; }
            dst(cc) += m_coefs[ii]*src(cc);
        }
    }
}

unsigned int MBPointInterpOp::size() const
{
    return m_coefs.size();
}

void MBPointInterpOp::print() const
{
    pr_out() << "--------------------------------" << std::endl;
    pr_out() << m_dst << std::endl;
    for (int ii = 0; ii < size(); ii++) {
        pr_out() << "\tcoef: " << m_coefs[ii] << " | " << m_srcs[ii] << std::endl;
    }
}

MBInterpOp::MBInterpOp()
{
    PR_TIME("MBInterpOp::constructor");
}

MBInterpOp::MBInterpOp(Point a_ghost, unsigned int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    m_order = a_order;
    m_ghost = a_ghost;
}

template<template<MemType> typename MAP, MemType MEM>
MBInterpOp::MBInterpOp(
        const MBLevelMap<MAP,MEM>&    a_map,
        Point a_ghost,
        int a_order)
{
    PR_TIME("MBInterpOp::constructor");
    define(a_map, a_ghost, a_order);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
                        const MBLevelMap<MAP,MEM>&    a_map,
                        Point a_ghost,
                        int a_order)
{
  PR_TIME("MBInterpOp::define1");
    MBInterpLayout interpLayout(a_map.layout(), a_order);
    define(a_map, interpLayout, a_ghost, a_order);
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
                        const MBLevelMap<MAP,MEM>&  a_map,
                        MBInterpLayout              a_interpLayout,
                        Point                       a_ghost,
                        int                         a_order)
{
    PR_TIME("MBInterpOp::define2");
    PROTO_ASSERT(a_interpLayout.layout().compatible(a_map.layout()),
        "MBInterpOp::define | Error: Incompatible layouts among inputs");

    m_order = a_order;
    m_ghost = a_ghost;
    m_interpLayout = a_interpLayout;
    for (auto iter : m_interpLayout.layout())
    {
        for (auto dir : Box::Kernel(1))
        {
           addPointOps(iter, dir, a_map);
        }
    }
}
template<MemType MEM>
void MBInterpOp::mapToPatch(MBDisjointBoxLayout a_layout)
{
  PR_TIME("MBInterpOp::mapToPatch");
  // PC: new version of apply using BoxData representation of stencils.
  int lsize = a_layout.localSize();
  m_srcSizes.resize(lsize);
  m_srcCoefs.resize(lsize);
  m_srcPoints.resize(lsize);
  m_srcDataPoints.resize(lsize);
  // cout << m_ghost << endl;
  for (auto dit : a_layout)
    {
      int localindex = dit.local();
      Box bx=a_layout[dit].grow(m_ghost);
      unsigned int block = a_layout.block(dit);
      Box bxblock = a_layout.domain().getBlock(block).box();
      Box bxvalid = a_layout[dit].grow(m_ghost);
      Point lowvalid = bxvalid.low();
      Point highvalid = bxvalid.high();
      lowvalid[0] = max(lowvalid[0],0);
      highvalid[0] = min(highvalid[0],bxblock.high()[0]);
      bxvalid = Box(lowvalid,highvalid);
      m_srcPoints[localindex] = make_shared<BoxData<Point,13,MEM> >(bxvalid);
      m_srcCoefs[localindex] = make_shared<BoxData<double,13,MEM> >(bxvalid);
      m_srcSizes[localindex] = make_shared<BoxData<uint,1,MEM> >(bxvalid);
      m_srcDataPoints[localindex] =
        make_shared<BoxData<MBDataPoint,13,MEM> >(bxvalid);
      
      m_srcCoefs[localindex]->setVal(0.);
      m_srcSizes[localindex]->setVal(0);
      m_srcPoints[localindex]->setVal(Point::Zeros());
    }
  // store stencils into BoxData representation.
  for (auto& op : m_ops)
    {
      auto mbdpt = op->target();
      auto index = mbdpt.index();
      
      auto pt = mbdpt.point();
      auto localindex = index.local();
      MBDataPoint mbdptsrc;
      Box bx = (*(m_srcSizes[localindex])).box();
     
      if (bx.low()[0] == pt[0])
        {
          
          auto thissize = op->coefs().size();
          (*(m_srcSizes[localindex]))(pt) = thissize;
          for (int cc = 0; cc < thissize; cc++)
            {
              (*(m_srcPoints[localindex]))(pt,cc) = op->sources()[cc].point();
              (*(m_srcDataPoints[localindex]))(pt,cc) = op->sources()[cc];
              auto coef = op->coefs()[cc];
              (*(m_srcCoefs[localindex]))(pt,cc) = coef;
            }
          
        }   
    }
  for (int localindex = 0; localindex < a_layout.localSize(); localindex++)
    {
      Operator::spreadSlice(*m_srcPoints[localindex]);
      Operator::spreadSlice(*m_srcCoefs[localindex]);
      Operator::spreadSlice(*m_srcSizes[localindex]);
      Operator::spreadSlice(*m_srcDataPoints[localindex]);
    }
  HDF5Handler h5;
  for (auto dit : a_layout)
    {
      PR_TIME("MBInterpOp::mapToPatch:output");
      
      uint globalindex = dit.global();
      uint localindex = dit.local();
      if (procID() == 0)
        {
          cout << "in iterator loop: "
             << globalindex << ", "
             << localindex << endl;
        }
      Box bxout = (*m_srcPoints[localindex]).box().grow(Point(0,-2,-2));
      if (procID() == 0)
        {
          Box bx=a_layout[dit].grow(m_ghost);
          cout << "box, box out = " <<bx << ", " << bxout << endl;
        }
      BoxData<double,26> points_bd(bxout);
      BoxData<double,1> sizes_bd(bxout);
      BoxData<double,13> coefs_bd(bxout);
      points_bd.setVal(0.);
      sizes_bd.setVal(0.);
      coefs_bd.setVal(0.);
      BoxData<Point,13>& ptsrc =*(m_srcPoints[localindex]);
      int count = 0;
      for (auto bit : bxout)
        {
          auto thissize = (*m_srcSizes[localindex])(bit);
          sizes_bd(bit) = thissize;
          for (int cc = 0; cc < thissize; cc++)
            {
              Point pt = ptsrc(bit,cc);
              points_bd(bit,cc) = pt[1];
              points_bd(bit,cc+13) = pt[2];
              coefs_bd(bit,cc) = (*m_srcCoefs[localindex])(bit,cc);            
            }
          count++;
        }
      auto block = a_layout.block(dit);
      
      // forallInPlace([] PROTO_LAMBDA(Var<double,26>& a_pointstoints,
      //                               const Var<Point,13>& a_ptsrcs)
      //               {
      //                 for (int cc = 0; cc < 13; cc++)
      //                   {
      //                     a_pointstoints(cc) = a_ptsrcs(cc)[1];
      //                     a_pointstoints(cc+13) = a_ptsrcs(cc)[2];
      //                   }
      //               },pointsToInts,ptsrc);
      //h5.writePatch(1.0,points_bd,"ptsources" + to_string(globalindex));
      //h5.writePatch(1.0,coefs_bd,"coefs" + to_string(globalindex));
      //h5.writePatch(1.0,sizes_bd,"sizes" + to_string(globalindex));
    }
}
template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::define(
        const std::vector<std::shared_ptr<MBLevelMap<MAP,MEM>>>&  a_maps,
        MBInterpLayout              a_interpLayout,
        Point                       a_ghost,
        int                         a_order)
{
    PR_TIME("MBInterpOp::define3");
    //cout << "define3 ghost input = " << a_ghost << endl;
    PROTO_ASSERT(a_maps.size() == a_interpLayout.layout().numBlocks(),
        "MBInterpOp::define | Error: Incorrect number of maps. There must be one map per block");
    for (auto bi = 0; bi < a_interpLayout.layout().numBlocks(); bi++)
    {
        PROTO_ASSERT(a_interpLayout.layout().compatible(a_maps[bi]->layout()),
            "MBInterpOp::define | Error: Incompatible layouts among inputs");
        PROTO_ASSERT(a_maps[bi]->ghost() == a_maps[0]->ghost(),
            "MBInterpOp::define | Error: Input maps have different ghost regions");
    }

    m_order = a_order;
    m_ghost = a_ghost;
    m_interpLayout = a_interpLayout;
    for (auto iter : m_interpLayout.layout())
    {
        auto block = m_interpLayout.layout().block(iter);
        auto map = a_maps[block];
        for (auto dir : Box::Kernel(1))
        {
           addPointOps(iter, dir, *map);
        }
    }
    // New storage scheme for stencils (PC).
    //mapToPatch<MEM>(m_interpLayout.layout());  
}

template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::addPointOps(const MBIndex& index, Point dir, const MBLevelMap<MAP,MEM>& map)
{
  PR_TIME("MBInterpOp::addPointOps");
  auto& layout = m_interpLayout.layout();
  PatchID patch = layout.point(index);
  PatchID adjPatch = patch + dir;
  BlockIndex block = layout.block(index);
  auto& patchDomain = layout.patchDomain();
  if (index.global() == 0)
    {
      // cout << "patch = " << patch
      //      << ", adjPatch = " << adjPatch
      //      << ", block = " << block << "\n"
      //      << "patchDomain.isPointInInterior(adjPatch, block) = "
      //      << patchDomain.isPointInInterior(adjPatch, block) << "\n"
      //      << "patchDomain.isPointInDomainBoundary(adjPatch, block) = " 
      //      << patchDomain.isPointInDomainBoundary(adjPatch, block)
      //      << endl;
    }
  if (patchDomain.isPointInInterior(adjPatch, block)) { return; }
  else if (patchDomain.isPointInDomainBoundary(adjPatch, block)) { return; }
  else {
    Box patchBox = layout[index];
    Point offset(0,-2,-2);
    Box boundBox = patchBox.adjacent((m_ghost + offset) * dir);
    // if (block == 0)
    //   {
    //     cout << endl;
    //     cout << "patchBox = " << patchBox << "\n"
    //          << "boundBox = " << boundBox << "\n"
    //          << "patch = " << patch << "\n"
    //          << "adjPatch = " << adjPatch << "\n"
    //          << "m_ghost = " << m_ghost 
    //          << endl;
    //   }
    for (auto boundPoint : boundBox)
      {
        if (this->containsOp(boundPoint, index))
          {
            continue;
          }
        addPointOp(boundPoint, index, map);
        // if (index.global() == 17)
        //   {
        //     cout << "boundPoint = " << boundPoint
        //          << ", patchBox = " <<  patchBox
        //          << ", boundBox = " << boundBox
        //          << endl;
        //   }
      }
  }
  
      
}

bool MBInterpOp::containsOp(Point dataPoint, const MBIndex& index) const
{
    PR_TIME("MBInterpOp::containsOp");
    MBPointID key(dataPoint, index.global());
    return (m_points.count(key) > 0);
}
template<template<MemType> typename MAP, MemType MEM>
void MBInterpOp::addPointOp(
                            Point dataPoint,
                            const MBIndex& index,
                            const MBLevelMap<MAP,MEM>& map)
{
    PR_TIME("MBInterpOp::addPointOp");
    PROTO_ASSERT(!containsOp(dataPoint, index),
        "MBInterpOp::addPointOp | Error: Trying to add duplicate op");
    auto layout = map.layout();
    //cout << "partition to block: " << layout.partitionPtr()->block(index) << endl;
    auto domain = layout.domain();
    auto destBlock = layout.block(index);
    auto destBlockBox = domain.getBlock(destBlock).box();
    if (destBlockBox.containsPoint(dataPoint)) return;
    MBDataPoint dstDataPoint(index, dataPoint, layout);
    // cout << "constructor: dstDataPoint  = " << dstDataPoint.point() <<
    //  " , " << dstDataPoint.block() << endl;    
    auto op =
      std::make_shared<MBPointInterpOp>(dstDataPoint, m_ghost, map, m_order);
    
    m_ops.push_back(op);
    MBPointID key(dataPoint, index.global());
    m_points.insert(key);

    copyAlongSymmetryAxes(*op);
}

void MBInterpOp::copyAlongSymmetryAxes(const MBPointInterpOp& pointOp)
{
    PR_TIME("MBInterpOp::copyAlongSymmetryAxis");
    Point dataPoint = pointOp.target().point();
    MBIndex index = pointOp.target().index();
    Box patchBox = m_interpLayout.layout()[index];
    Point boundaryDir = patchBox.whichBoundaryContains(dataPoint);
    Box patchBoundBox = patchBox.adjacent(boundaryDir * m_ghost);
    for (int dd = 0; dd < DIM; dd++)
    {
        PR_TIME("MBInterpOp::copyAlongSymmetryAxis");

        if (!m_interpLayout.isCopied(dd)) { continue; }

        Point low = dataPoint;
        Point high = dataPoint;
        low[dd] = patchBoundBox.low()[dd];
        high[dd] = patchBoundBox.high()[dd];
        Box axisBox(low, high);
        PROTO_ASSERT(axisBox.containsPoint(dataPoint),
            "MBInterpOp::copyAlongSymmetryAxes | Error: Data Corruption");
        for (auto axisPoint : axisBox)
        {
            if (containsOp(axisPoint, index)) { continue; }
            Point shift = axisPoint - dataPoint;
            addShiftedOp(pointOp, shift);
        }
    }
}
void MBInterpOp::addShiftedOp(
    const MBPointInterpOp& baseOp,
    Point shift)
{
    PR_TIME("MBInterpOp::addShiftedOp");
    PROTO_ASSERT(containsOp(baseOp.target().point(), baseOp.target().index()),
        "MBInterpOp::adPointOp | Error: Trying to add duplicate op");
    if (shift == Point::Zeros()) { return; }
    MBDataPoint shiftDst = baseOp.target().shift(shift);
    std::vector<MBDataPoint> shiftSrc;
    for (auto si : baseOp.sources())
    {
        shiftSrc.push_back(si.shift(shift));
    }
    auto coefs = baseOp.coefs();
    auto shiftedOp = std::make_shared<MBPointInterpOp>(shiftDst, shiftSrc, coefs);
    
    m_ops.push_back(shiftedOp);
    MBPointID key(shiftDst.point(), baseOp.target().index().global());
    m_points.insert(key);
}

void MBInterpOp::insert(std::shared_ptr<MBPointInterpOp> a_op)
{
    PR_TIME("MBInterpOp::insert");
    Point dst = a_op->target().point();
    int globalIndex = a_op->target().index();
    MBPointID key(dst, globalIndex);
    PROTO_ASSERT(m_points.count(key) == 0,
            "MBInterpOp::insert | Tried to insert a duplicate MBPointInterpOp");

    m_ops.push_back(a_op);
    m_points.insert(key);
}

template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::applyNew(
                       MBLevelBoxData<T, C, MEM>& a_dst,
                       MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBInterpOp::applyNew");
    
    auto layout = a_dst.layout();
    for (auto dit : layout)
      {
        auto indexlocal = dit.local();
        unsigned int block = layout.block(dit);
        Box bxblock = layout.domain().getBlock(block).box();
        Box bxapply = a_dst[dit].box();
        auto high = bxapply.high();
        auto low = bxapply.low();
        high[0] = min(high[0],bxblock.high()[0]);
        low[0] = max(low[0],bxblock.low()[0]);
        bxapply = Box(low,high);
        for (auto bit : bxapply)
          {
            auto size = (*m_srcSizes[dit])(bit);
            if ((size > 0) && !bxblock.containsPoint(bit))
              {
                auto dest = a_dst[dit].var(bit);
                auto mbdpts = m_srcDataPoints[indexlocal]->var(bit);
                auto coefs = m_srcCoefs[indexlocal]->var(bit);
                for (uint ii = 0; ii < size; ii++)
                  {
                    auto sources = a_src[mbdpts(ii)];
                    if (ii == 0)
                      {
                        for (uint cc = 0; cc < C; cc++)
                          {
                            dest(cc) = 0.;
                          }
                      }
                    for (uint cc = 0; cc < C; cc++)
                      {
                        dest(cc) += coefs(ii)*sources(cc);
                      }
                  }
              }
          }
      }
    // if ((!a_src[dit].box().containsPoint(points(ii)))
    //     && (block == 0))
    //   {
    //     cout << "src box = " << a_src[dit].box() << "\n"
    //          << "dest box = " << a_dst[dit].box() << "\n"
    //          << "bxapply = " << bxapply << "\n"
    //          << "dest point = " << bit << "\n"
    //          << "src point = " << points(ii) << "\n"
    //          << "coef # = " << ii << "\n"
    //          << endl << endl;
    // for (uint iidebug = 0; iidebug < size; iidebug++)
    //     {
    //       cout << "ii = " << iidebug 
    //       << ", src point (ii) = " << points(iidebug) 
    //       << ", coef (ii) = " << coefs(iidebug)
    //        << endl;
    //     }
    //}
}
    template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::apply(
                       MBLevelBoxData<T, C, MEM>& a_dst,
                       MBLevelBoxData<T, C, MEM>& a_src)
{
    PR_TIME("MBInterpOp::apply");
    for (auto& op : m_ops)
    {      
      // Apply interpolation operator only where a_dst is defined.
        Point pt = op->target().point();
        auto index = op->target().index();
        if (a_dst[index].box().containsPoint(pt))
          op->apply(a_dst, a_src);
    }
}
void MBInterpOp::print() const
{
    pr_out() << "Number of PointInterpOps: " << m_ops.size() << std::endl;
    for (auto& op : m_ops)
    {
        op->print(); 
    }
}
// void MBInterpOp::copyAxis(unsigned int a_axis, bool a_flag)
// {
//     m_copyAxes[a_axis] = a_flag;
// }
template<typename T, unsigned int C, MemType MEM>
void MBInterpOp::printErrorPoints(
        MBLevelBoxData<T, C, MEM>& a_err,
        T                          a_eps)
{
    auto& layout = a_err.layout();
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& err_i = a_err[iter];
        for (auto pi : err_i.box())
        {
            for (int cc = 0; cc < C; cc++)
            {
                if (err_i(pi, cc) > a_eps)
                {   
                    for (auto& op : m_ops)
                    {
                        auto& dst = op->target();
                        if (dst.point() == pi && dst.dstBlock() == block)
                        {
                            std::cout << "\nError Point found | error: " << err_i(pi, cc) << " | target: " << pi << " | block: " << block << std::endl;
                            std::cout << "Source Points: " << std::endl;
                            for (auto& si : op->sources())
                            {
                                std::cout << "\t" << si << std::endl;
                            }
                            break;
                        }
                    }
                    break;
                }
                
            }
        }
    }
}


template<template<MemType> typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(MBLevelBoxData<T, C, MEM, CTR>& a_data, unsigned int a_order)
{
    MBLevelMap<MAP, MEM> map(a_data.layout(), a_data.ghost());
    interpBoundaries<MAP, T, C, MEM, CTR>(a_data, map, a_order);
}

template<template<MemType> typename MAP, typename T, unsigned int C, MemType MEM, Centering CTR>
void interpBoundaries(
        MBLevelBoxData<T, C, MEM, CTR>& a_data,
        MBLevelMap<MAP, MEM>& a_map,
        unsigned int a_order)
{
    std::vector<Point> footprint;
    Box K = Box::Kernel(2);
    for (auto pi : K)
    {
        if (pi.abs().sum() <= 2)
        {
            footprint.push_back(pi);
        }
    }
    MBInterpOp op(a_data.ghost()[0], a_order);
    op.define(a_map, footprint);
    op.apply(a_data, a_data);
}
