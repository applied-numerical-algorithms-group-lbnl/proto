
template<MemType MEM>
MBMapOp<MEM>::MBMapOp(Array<double, DIM> a_dx, unsigned int a_block)
{
    define(a_dx, a_block);
}

template<MemType MEM>
void MBMapOp<MEM>::define(Array<double, DIM> a_dx, unsigned int a_block)
{
    m_dx = a_dx;
    m_block = a_block;
    m_c2c = Stencil<double>::CornersToCells(4);
    this->init();
}

template<MemType MEM>
void MBMapOp<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        unsigned int a_block)
{
    FluxBoxData<double, DIM, MEM> NT(a_J.box());
    apply(a_X, a_J, NT, a_block);
}

template<MemType MEM>
void MBMapOp<MEM>::apply(
        BoxData<double, DIM, MEM>& a_X,
        BoxData<double, 1, MEM>& a_J,
        FluxBoxData<double, DIM, MEM>& a_NT,
        unsigned int a_block)
{
    apply(a_X, a_J, a_block);
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
