
MBMapOp::MBMapOp(Array<double, DIM> a_dx, unsigned int a_block = 0)
{
    define(a_dx, a_block);
}

void MBMapOp::define(Array<double, DIM> a_dx, unsigned int a_block = 0)
{
    m_dx = a_dx;
    m_block = a_block;
    m_c2c = Stencil<double>::CornersToCells(4);
    this->init();
}

void MBMapOp::init() {}

template<MemType MEM>
void MBMapOp::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_block)
{
    FluxBoxData<double, DIM, MEM> NT(a_J.box());
    apply(a_X, a_J, NT, a_block);
}

template<MemType MEM>
void MBMapOp::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT,
        unsigned int a_block)
{
    apply(a_X, a_J, a_block);
}


template<MemType MEM>
void MBMapOp::doApply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_computeBlock,
        unsigned int a_outBlock)
{
    if (a_computeBlock != a_outBlock)
    {
        const auto& domain = layout().domain();
        Box BX = domain.convert(a_X.box(), a_outBlock, a_computeBlock, PR_NODE);
        Box BJ = domain.convert(a_J.box(), a_outBlock, a_computeBlock, PR_CELL);
        BoxData<double, DIM, MEM> X_tmp(BX);
        BoxData<double, 1, MEM> J_tmp(BJ);
        apply(X_tmp, J_tmp, a_computeBlock);
        Point conn = domain.graph().connectivity(a_computeBlock, a_outBlock);
        auto R = domain.graph().rotation(a_computeBlock, conn, a_outBlock);
        X_tmp.copyTo(a_X, R);
        J_tmp.copyTo(a_J, R);
    } else {
        apply(a_X, a_J, a_computeBlock);
    }
}

template<MemType MEM>
BoxData<double, DIM, MEM> MBMapOp::cellCentered(
        const Box& a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock)
{
    auto X0 = cellAveraged(a_box.grow(1), a_computeBlock, a_outBlock);
    return Operator::deconvolve(X0);
}

template<MemType MEM>
BoxData<double, DIM, MEM> MBMapOp::cellAveraged(
        const Box& a_box,
        unsigned int a_computeBlock,
        unsigned int a_outBlock)
{
    Box B0 = m_c2c.domain(a_box);
    BoxData<double, DIM, MEM> X0(B0);
    BoxData<double, 1, MEM> J0(B0.extrude(Point::Ones(), -1)); //unused
    doApply(X0, J0, a_computeBlock, a_outBlock);
    return m_c2c(X0);
}

Array<double, DIM> MBMapOp::cellCentered(const MBDataPoint& a_point)
{
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.point, a_point.point);
    auto X = cellCentered(B0, computeBlock, domainBlock);
    return X.array(a_point.point);
}

Array<double, DIM> MBMapOp::cellAveraged(const MBDataPoint& a_point)
{
    auto computeBlock = a_point.srcBlock();
    auto domainBlock = a_point.dstBlock();
    Box B0(a_point.point, a_point.point);
    auto X = cellAveraged(B0, computeBlock, domainBlock);
    return X.array(a_point.point);
}

const Array<double, DIM>& MBMapOp::dx() const { return m_dx; }
}

unsigned int MBMapOp::block() const { return m_block; }

template<MemType MEM>
BoxData<double, DIM, MEM> MBMapOp::X(const Box& a_box, const Array<double, DIM>& a_dx) const
{
    BoxData<double, DIM, MEM> XOut(a_box);
    forallInPlace_p(
        [] PROTO_LAMBDA (Point& a_pt, Var<double, DIM, MEM>& a_x, Array<double, DIM> a_dx)
        {
            for (int dir = 0; dir < DIM; dir++)
            {
                a_x(dir) = a_pt[dir]*a_dx[dir];
            }
        }, XOut, a_dx);
    return XOut;
}
