
template<MemType MEM>
MBMapOp<MEM>::MBMapOp(
        Array<double, DIM> a_dx,
        unsigned int a_block,
        unsigned int a_focalBlock)
{
    define(a_dx, a_block, a_focalBlock);
}

template<MemType MEM>
void MBMapOp<MEM>::define(
        Array<double, DIM> a_dx,
        unsigned int a_block,
        unsigned int a_focalBlock)
{
    m_dx = a_dx;
    m_block = a_block;
    m_focalBlock = a_focalBlock;
    m_c2c = Stencil<double>::CornersToCells(4);
    this->init();
}

//NB: these two functions don't form a loop because the user must define one of them
// in the derived class, hence breaking the loop and allowing the undefined version
// to funnel into the defined one. 
template<MemType MEM>
void MBMapOp<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J)
{
    FluxBoxData<double, DIM, MEM> NT(a_J.box());
    apply(a_X, a_J, NT);
}

//NB: these two functions don't form a loop because the user must define one of them
// in the derived class, hence breaking the loop and allowing the undefined version
// to funnel into the defined one. 
template<MemType MEM>
void MBMapOp<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT)
{
    apply(a_X, a_J);
}

template<MemType MEM>
const Array<double, DIM>& MBMapOp<MEM>::dx() const { return m_dx; }

template<MemType MEM>
unsigned int MBMapOp<MEM>::block() const { return m_block; }

template<MemType MEM>
BoxData<double, DIM, MEM> MBMapOp<MEM>::X(const Box& a_box, const Array<double, DIM>& a_dx) const
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
