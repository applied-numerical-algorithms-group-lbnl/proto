//==========================================================================================
// CONVOLVE
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    ConvolveOp<T> op;
    op(a_avg, a_ctr, a_2nd);
    // PR_TIME("Operator::convolve");

    // auto L = Stencil<T>::Laplacian();
    // Stencil<T> I = 1.0*Shift::Zeros();
    // Box range = L.range(a_2nd.box());
    
    // a_avg |= L(a_2nd, range, 1.0/24.0);
    // a_avg += I(a_ctr, range);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_convolve(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().containsBox(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,MEM,D,E> a_avg(a_ctr.box() & a_2nd.box().grow(-1));
    _convolve(a_avg, a_ctr, a_2nd);
    return a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
convolve(const BoxData<T,C,MEM,D,E>& a_ctr)
{
    return _convolve(a_ctr, a_ctr);
}

template<class T, unsigned int C, MemType MEM, Centering CTR>
void
_convolve(
        LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_2nd)
{
    //TODO: switch functionality based on value of CTR
    PROTO_ASSERT(CTR == PR_CELL,
            "Operator::_convolve | Error: Operator not implemented for this centering.");
    for (auto iter = a_avg.begin(); iter.ok(); ++iter)
    {
        auto& avg =  a_avg[*iter];
        auto& ctr =  a_ctr[*iter];
        auto& ctr2 = a_2nd[*iter];
        _convolve(avg, ctr, ctr2);
    }
}

template<class T, unsigned int C, MemType MEM, Centering CTR>
void
convolve(
        LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}

template<class T, unsigned int C, MemType MEM, Centering CTR,
    typename Func, typename... Srcs>
void
initConvolve(
        LevelBoxData<T,C,MEM,CTR>& a_data,
        Func& a_func,
        Srcs... a_srcs)
{
    for (auto iter : a_data.layout())
    {
        auto& data = a_data[iter];
        BoxData<T,C,MEM> tmp(data.box().grow(1));
        forallInPlace_p(a_func, tmp, a_srcs...);
        Operator::convolve(data, tmp);
    }
}

#ifdef PR_AMR
template<class T, unsigned int C, MemType MEM, Centering CTR,
    typename Func, typename... Srcs>
void
initConvolve(
        AMRData<T,C,MEM,CTR>& a_data,
        double a_dx0,
        Func& a_func,
        Srcs... a_srcs)
{
    //TODO: Make this anisotropic
    double dx = a_dx0;
    for (int lvl = 0; lvl < a_data.numLevels(); lvl++)
    {
        Operator::initConvolve(a_data[lvl], a_func, dx, a_srcs...);
        if (lvl < a_data.numLevels()-1)
        {
            dx /= a_data.grid().refRatio(lvl)[0];
        }
    }
}
#endif

//==========================================================================================
// DECONVOLVE
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    // PR_TIME("Operator::deconvolve");
    // PROTO_ASSERT(a_avg.box().containsBox(a_2nd.box().grow(-1)),
    //         "error Operator::convolve | centered data defined on too small a box.");

    // a_ctr |= Stencil<T>::Laplacian()(a_2nd, -1.0/24);
    // a_ctr += a_avg;
    DeconvolveOp<T> op;
    op(a_ctr, a_avg, a_2nd);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_deconvolve(
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().containsBox(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,MEM,D,E> a_ctr(a_avg.box() & a_2nd.box().grow(-1));
    _deconvolve(a_ctr, a_avg, a_2nd);
    return a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
deconvolve(const BoxData<T,C,MEM,D,E>& a_avg)
{
    return _deconvolve(a_avg, a_avg);
}

template<class T, unsigned int C, MemType MEM, Centering CTR>
void
_deconvolve(
        LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_2nd)
{
    //TODO: switch functionality based on value of CTR
    for (auto iter = a_avg.begin(); iter.ok(); ++iter)
    {
        auto& ctr =  a_ctr[*iter];
        auto& avg =  a_avg[*iter];
        auto& ctr2 = a_2nd[*iter];
        deconvolve(ctr, avg, ctr2);
    }
}

template<class T, unsigned int C, MemType MEM, Centering CTR>
void
deconvolve(
        LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}

//==========================================================================================
// CONVOLVE FACE
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    ConvolveFaceOp<T> op(a_dir);
    op(a_avg, a_ctr, a_2nd);
    // PR_TIME("Operator::convolveFace");
    // PROTO_ASSERT(a_ctr.box().containsBox(a_2nd.box().grow(-1).grow(a_dir,1)),
    // "Error in Operator::convolveFace | Insufficient source data.");
    // a_avg |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, 1.0/24.0);
    // a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_convolveFace(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    auto b = a_2nd.box().grow(-1).grow(a_dir, 1);
    b &= a_ctr.box();
    BoxData<T,C,MEM,D,E> avg(b);
    _convolveFace(avg, a_ctr, a_2nd, a_dir);
    return avg;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    _convolveFace(a_avg, a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
convolveFace(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    return _convolveFace(a_ctr, a_ctr, a_dir);
}

//==========================================================================================
// DECONVOLVE FACE
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    // PR_TIME("Operator::deconvolveFace");
    // PROTO_ASSERT((a_2nd.box().grow(-1).grow(a_dir,1)).containsBox(a_ctr.box()) &&
    //              a_avg.box().containsBox(a_ctr.box()),
    // "Error in Operator::convolveFace | Insufficient source data.");
    // a_ctr |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, -1.0/24.0);
    // a_ctr += a_avg;
    DeconvolveFaceOp<T> op(a_dir);
    op(a_ctr, a_avg, a_2nd);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_deconvolveFace(
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    auto b = a_2nd.box().grow(-1).grow(a_dir, 1);
    b &= a_avg.box();
    BoxData<T,C,MEM,D,E> ctr(b);
    _deconvolveFace(ctr, a_avg, a_2nd, a_dir);
    return ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    _deconvolveFace(a_ctr, a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
deconvolveFace(
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    return _deconvolveFace(a_avg, a_avg, a_dir);
}

//==========================================================================================
// CELL PRODUCT
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_cellProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    auto b = (a_p2nd.box().grow(-1) & a_q2nd.box().grow(-1));
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,MEM,D,E> pq(b);
    _cellProduct(pq, a_pavg, a_qavg, a_p2nd, a_q2nd);
    return pq;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    PR_TIME("Operator::cellProduct");
    CellAverageProductOp<T> op;
    op(a_product, a_pavg, a_qavg, a_p2nd, a_q2nd);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
cellProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellProduct(a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellProduct(a_product, a_pavg, a_qavg, a_pavg, a_qavg);
}

//==========================================================================================
// FACE PRODUCT
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_faceProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    Proto::Box b = (a_p2nd.box() & a_q2nd.box()).grow(-1).grow(a_dir, 1);
    
    BoxData<T,C,MEM,D,E> out(b);
    _faceProduct(out, a_pavg, a_qavg, a_p2nd, a_q2nd, a_dir);
    return out;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    FaceAverageProductOp<T> op(a_dir);
    op(a_dest, a_pavg, a_qavg, a_p2nd, a_q2nd);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
faceProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceProduct(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceProduct(a_dest, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

//==========================================================================================
// CELL QUOTIENT
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_cellQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    auto b = (a_p2nd.box().grow(-1) & a_q2nd.box().grow(-1));
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,MEM,D,E> quotient(b);
    _cellQuotient(quotient, a_pavg, a_qavg, a_p2nd, a_q2nd);
    return quotient;

}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    PR_TIME("Operator::cellQuotient");
    PROTO_ASSERT(a_pavg.box().containsBox(a_p2nd.box().grow(-1)),
            "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().containsBox(a_q2nd.box().grow(-1)),
            "Operator::cellQuotient | incompatible patch sizes.");
    CellAverageQuotientOp<T> op;
    op(a_quotient, a_pavg, a_qavg, a_p2nd, a_q2nd);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
cellQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellQuotient(a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg);
}

//==========================================================================================
// FACE QUOTIENT
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_faceQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    Proto::Box b = (a_p2nd.box() & a_q2nd.box()).grow(-1).grow(a_dir, 1);
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,MEM,D,E> out(b);
    _faceQuotient(out, a_pavg, a_qavg, a_p2nd, a_q2nd, a_dir);
    return out;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
_faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    PR_TIME("Operator::faceQuotient");
    PROTO_ASSERT(a_pavg.box().containsBox(a_p2nd.box().grow(-1)),
            "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().containsBox(a_q2nd.box().grow(-1)),
            "Operator::cellQuotient | incompatible patch sizes.");
    FaceAverageQuotientOp<T> op(a_dir);
    op(a_quotient, a_pavg, a_qavg, a_p2nd, a_q2nd);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
faceQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceQuotient(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void
faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

//==========================================================================================
// FACE TENSOR PRODUCT

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
BoxData<T,C,MEM,D,E>
_faceTensorProduct(
        const BoxData<T,C,MEM,D,E>& a_ql4,
        const BoxData<T,1,MEM,1,1>& a_qr4,
        const BoxData<T,C,MEM,D,E>& a_ql2,
        const BoxData<T,1,MEM,1,1>& a_qr2,        
        int a_dir)
{
    Box bx = a_ql4.box()&a_qr4.box();
    for (int dir = 0; dir < DIM; dir++)
    {
        if (dir!=a_dir)
        {
            bx = bx & a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
        }
    }
    BoxData<T,C,MEM,D,E> quot(bx);
    for (int ee = 0; ee < E;ee++)
    for (int dd = 0; dd < D; dd++)
    for (int cc = 0; cc < C; cc++)
    {
        auto quotcomp = slice(quot,cc,dd,ee);
        auto ql4comp = slice(a_ql4,cc,dd,ee);
        auto ql2comp = slice(a_ql2,cc,dd,ee);
        Operator::_faceProduct(quotcomp,ql4comp,a_qr4,ql2comp,a_qr2,a_dir);
    }
    return quot;
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
void
_faceTensorProduct(
        BoxData<T,C,MEM,D,E>&       a_product,
        const BoxData<T,C,MEM,D,E>& a_ql4,
        const BoxData<T,1,MEM,1,1>& a_qr4,
        const BoxData<T,C,MEM,D,E>& a_ql2,
        const BoxData<T,1,MEM,1,1>& a_qr2,
        int a_dir)
{
    MayDay<void>::Abort("Operator::_faceTensorProduct | Error: Not Implemented.");
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
BoxData<T,C,MEM,D,E>
faceTensorProduct(
        const BoxData<T,C,MEM,D,E>& a_ql,
        const BoxData<T,1,MEM,1,1>& a_qr,
        int a_dir)
{
    MayDay<void>::Abort("Operator::faceTensorProduct | Error: Not Implemented.");
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
void
faceTensorProduct(
        BoxData<T,C,MEM,D,E>&       a_product,
        const BoxData<T,C,MEM,D,E>& a_ql,
        const BoxData<T,1,MEM,1,1>& a_qr,
        int a_dir)
{
    MayDay<void>::Abort("Operator::faceTensorProduct | Error: Not Implemented.");
}

//==========================================================================================
// FACE TENSOR QUOTIENT

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
BoxData<T,C,MEM,D,E>
_faceTensorQuotient(
        const BoxData<T,C,MEM,D,E>& a_ql4,
        const BoxData<T,1,MEM,1,1>& a_qr4,
        const BoxData<T,C,MEM,D,E>& a_ql2,
        const BoxData<T,1,MEM,1,1>& a_qr2,        
        int a_dir)
{
    PR_TIME("Operator::faceTensorQuotient");
    FaceAverageTensorQuotientOp<T> op(a_dir);
    return op(a_ql4, a_qr4, a_ql2, a_qr2);
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
void
_faceTensorQuotient(
        BoxData<T,C,MEM,D,E>&       a_quotient,
        const BoxData<T,C,MEM,D,E>& a_ql4,
        const BoxData<T,1,MEM,1,1>& a_qr4,
        const BoxData<T,C,MEM,D,E>& a_ql2,
        const BoxData<T,1,MEM,1,1>& a_qr2,
        int a_dir)
{
    MayDay<void>::Abort("Operator::_faceTensorQuotient | Error: Not Implemented.");
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
BoxData<T,C,MEM,D,E>
faceTensorQuotient(
        const BoxData<T,C,MEM,D,E>& a_ql,
        const BoxData<T,1,MEM,1,1>& a_qr,
        int a_dir)
{
    MayDay<void>::Abort("Operator::faceTensorQuotient | Error: Not Implemented.");
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
void
faceTensorQuotient(
        BoxData<T,C,MEM,D,E>&       a_quotient,
        const BoxData<T,C,MEM,D,E>& a_ql,
        const BoxData<T,1,MEM,1,1>& a_qr,
        int a_dir)
{
    MayDay<void>::Abort("Operator::faceTensorQuotient | Error: Not Implemented.");
}

//==========================================================================================
// COFACTOR

template <class T, MemType MEM>
void
cofactor(
        BoxData<T,DIM,MEM>&       a_N,
        const BoxData<T,DIM,MEM>& a_X,
        int a_dir)
{
    MayDay<void>::Abort("Operator::cofactor | Error: Not implemented.");
}
template <class T, MemType MEM>
BoxData<T,DIM,MEM>
cofactor(
        const BoxData<T,DIM,MEM>& a_X,
        int a_dir,
        Array<T,DIM> a_dx)
{
    auto N = cofactor(a_X, a_dir);
    T area = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        if (d == a_dir) {continue; }
        area *= a_dx[d];
    }
    N /= area;
    return N;
}

template <class T, MemType MEM>
BoxData<T,DIM,MEM>
cofactor(
        const BoxData<T,DIM,MEM>& a_X,
        int a_dir)
{
    // compute cofactor averaged on face dir from node centered X
    PR_TIME("Operator::cofactor");
#if DIM==2
    int dirperp = (a_dir+1)%DIM;
    Box bxret = a_X.box()&a_X.box().shift(dirperp,-1);
    BoxData<T,DIM,MEM> ret(bxret);
    ret.setToZero();

    for (int comp = 0;comp < DIM;comp++)
    {
        int compperp = (comp+1)%DIM;
        auto retslice = slice(ret,comp);
        auto Xslice = slice(a_X,compperp);
        T sign = 1.0;
        if (compperp!=dirperp) sign = -1.0;
        retslice += Stencil<T>::FluxDivergence(dirperp)(Xslice,sign);
    }
    return ret;
#elif DIM==3
    int dirperp1 = (a_dir+1)%DIM;.
    int dirperp2 = (a_dir+2)%DIM;

    // average X in direction 1
    // XAv1 is on edges parallel to dirperp1
    BoxData<T,DIM,MEM> XAv1 = Stencil<T>::faceToCell(dirperp1,4)(a_X);
    // dX/dx1
    // dXdxi1 is on edges parellel to dirperp1, 2nd order, undivided
    BoxData<T,DIM,MEM> dXdxi1 = Stencil<T>::FluxDivergence(dirperp1)(a_X);
    // dX/dx1 ^ X
    auto dXdxi1byX = Operator::_edgeCrossProduct3D(dXdxi1,XAv1,dXdxi1,XAv1,dirperp1);
    // average X in direction 2
    BoxData<T,DIM,MEM> XAv2 = Stencil<T>::faceToCell(dirperp2,4)(a_X);
    // dX/dx2
    BoxData<T,DIM,MEM> dXdxi2 = Stencil<T>::FluxDivergence(dirperp2)(a_X);
    // X ^ dX/dx2
    auto XbydXdxi2 = Operator::_edgeCrossProduct3D(XAv2,dXdxi2,XAv2,dXdxi2,dirperp2);
    // d/dx2( dX/dx1 ^ x)
    BoxData<T,DIM,MEM> d12 = Stencil<T>::FluxDivergence(dirperp2)(dXdxi1byX,.5);
    // d/dx1(X ^ dX/dx2)
    BoxData<T,DIM,MEM> d21 = Stencil<T>::FluxDivergence(dirperp1)(XbydXdxi2,.5);   
    // return  0.5*d/dx2(dX/dx1 ^ X) + 0.5*d/dx1(dX/dx2 ^ X)
    auto ret = forall<T,DIM,MEM>(
        []PROTO_LAMBDA(Var<T,DIM,MEM>& a_ret, Var<T,DIM,MEM>& a_d21,Var<T,DIM,MEM>& a_d12)
        {
            for (int comp = 0; comp < DIM; comp++)
            {
                a_ret(comp) = a_d12(comp) + a_d21(comp);
            }
        },d21,d12);
    
    return ret;
#else
    std::cout << "cofactor only implemented for DIM = 2,3" << endl;
    abort();
#endif
}

//==========================================================================================
// COFACTOR MATRIX

template<class T, MemType MEM>
void
cofactorMatrix(
        BoxData<T,DIM,MEM,DIM>& a_NTMatrix,
        const FluxBoxData<T,DIM,MEM>& a_NTFaces,
        //const Array<BoxData<T,DIM,MEM>, DIM>& a_NTFaces,
        int a_dirNorm)                         
{
    MayDay<void>::Abort("Operator::cofactorMatrix | Error: Function not implemented.");
}

// Copy lowest 0-index slice to all the other slices.
template<class T,unsigned int C,MemType MEM,unsigned int D, unsigned int E>
void spreadSlice(BoxData<T,C,MEM,D,E>& a_bd)
{
  PR_TIMERS("spreadSlice");
  Box bx= a_bd.box();
  Point low = bx.low();
  Point high = bx.high();
  Point low0 = low;
  Point high0 = high;
  high0[0] = low0[0];
  Box bx0(low0,high0);
  /* for (int c = 0; c < C; c++)
    {
      for (auto bit : bx0)
        {
          Point plow = bit;
          Point phigh = bit;
          plow[0]++;
          phigh[0] = high[0];
          Box pencil(plow,phigh);
          T baseval = *a_bd.data(bit,c);
          forallInPlace([]PROTO_LAMBDA
                        (Var<T,C,MEM>& a_data,
                         T a_baseval,
                         int a_c)
                        {
                          a_data(a_c) = a_baseval;
                        },pencil,a_bd,baseval,c);
        }
        }*/
                       
  // The version below is about 6x faster on a CPU than the forall version above.
  // Reduces the run time of CubedSphereGeometry by about 25%.
  for (int e = 0; e < E; e++)
  for (int d = 0; d < D; d++)
  for (int c = 0; c < C; c++)
    {
      for (auto bit : bx0)
        {
          T* data0 = a_bd.data(bit,c,d,e);
          for (int i = 1; i < bx.high()[0] - bx.low()[0] + 1; i++)
            {
              *(data0 + i) = *data0;
            }
        }
    }   
}

template<class T, MemType MEM>
void
cubedSphere_nSphere(
                    BoxData<T,DIM,MEM>& a_nSphere,
                    T a_h)                      
{
  // compute the location of the points in a_bx projected back to the surface of the
  // sphere of radius 1.
  PR_TIMERS("cubedSphere_nSphere");
  // HDF5Handler h5;
  Box bx = a_nSphere.box();
  Point low0 = bx.low();
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  // Fill in loweest slice using analytic formula for n.
  forallInPlace_p([]PROTO_LAMBDA
                 (Point a_pt,
                  Var<T,DIM,MEM>& a_nSphere,                 
                  T a_h)
                 {
                   T alpha = -.25*M_PI + .5*M_PI*a_pt[1]*a_h;
                   T beta = -.25*M_PI + .5*M_PI*a_pt[2]*a_h;
                   T tana = tan(alpha);
                   T tanb = tan(beta);
                   T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
                   a_nSphere(0) = norm;
                   a_nSphere(1) = tana*norm;
                   a_nSphere(2) = tanb*norm;
                 },bx0,a_nSphere,a_h); 
  spreadSlice(a_nSphere);
}                   
  // Compute the average of D_1(n) x n along 1-edges.
template<class T, MemType MEM>
void
cubedSphere_edge1(
                  BoxData<T,DIM,MEM>& a_d1n_cross_n,
                  T a_h)
{
  PR_TIMERS("cubedSphere_edge1");
#if QUADRATURE1D == 5
  Array<T,5> weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
#endif
#if QUADRATURE1D == 3
  Array<T,3> weights={1.0/6.0,2.0/3.0,1.0/6.0};
#endif
#if QUADRATURE1D == 2
  Array<T,2> weights={.5,.5};
#endif
  int q1d = QUADRATURE1D;
  Box bx = a_d1n_cross_n.box();
  Point low0 = bx.low();
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  forallInPlace_p([q1d]PROTO_LAMBDA
                  (Point a_pt,
                   Var<T,DIM,MEM>& a_d1n_cross_n,
                   Array<T,QUADRATURE1D>& a_weights,
                   T a_h)
                  {
                    a_d1n_cross_n(0) = 0.0;
                    a_d1n_cross_n(1) = 0.0;
                    a_d1n_cross_n(2) = 0.0;
                    T beta = -.25*M_PI + .5*M_PI*a_pt[2]*a_h ;
                    T tanb = tan(beta);
                    for (int k = 0; k < q1d;k++)
                      {
                        T alpha = -.25*M_PI + .5*M_PI*(a_pt[1]*a_h + k*a_h/(q1d-1));
                        T tana = tan(alpha);
                        T seca = 1.0/cos(alpha);
                        T normsq = 1.0/(tana*tana + tanb*tanb + 1.0);               
                        Array<T,DIM> cross;
                        cross[0] = seca*seca*tanb*normsq;
                        //cross[1] = 0;
                        cross[2] = -seca*seca*normsq;
                        a_d1n_cross_n(0) += cross[0]*a_weights[k];
                        //a_d1n_cross_n(1) += cross[1]*a_weights[k];
                        a_d1n_cross_n(2) += cross[2]*a_weights[k];
                      }
                  },bx0,a_d1n_cross_n,weights,a_h);
  // copy lowest slice to all the other slices. 
  spreadSlice(a_d1n_cross_n);
}
  // Compute the integral of n x D_2(n) along 2-edges.
template<class T, MemType MEM>
void
cubedSphere_edge2(
                  BoxData<T,DIM,MEM>& a_n_cross_d2n,
                  T a_h)
{
  PR_TIMERS("cubedSphere_edge2");
#if QUADRATURE1D == 5
  Array<T,5> weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
#endif
#if QUADRATURE1D == 3
  Array<T,3> weights={1.0/6.0,2.0/3.0,1.0/6.0};
#endif
#if QUADRATURE1D == 2
  Array<T,2> weights={.5,.5};
#endif
  int q1d = QUADRATURE1D;
  Box bx = a_n_cross_d2n.box();
  Point low0 = bx.low();
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  forallInPlace_p([q1d]PROTO_LAMBDA
                  (Point a_pt,
                   Var<T,DIM,MEM>&  a_n_cross_d2n,
                   Array<T,QUADRATURE1D>& weights,
                   T a_h)
                  {
                    a_n_cross_d2n(0) = 0.0;
                    a_n_cross_d2n(1) = 0.0;
                    a_n_cross_d2n(2) = 0.0;
                    T alpha = -.25*M_PI + .5*M_PI*a_pt[1]*a_h;
                    T tana = tan(alpha);
                    for (int k = 0; k < q1d;k++)
                      {
                        T beta = -.25*M_PI + .5*M_PI*(a_pt[2]*a_h + k*a_h/(q1d-1));
                        T tanb = tan(beta);
                        T secb =  1.0/cos(beta);
                        T normsq = 1.0/(tana*tana + tanb*tanb + 1.0);         
                        Array<T,DIM> cross;
                        cross[0] = tana*secb*secb*normsq;
                        cross[1] = -secb*secb*normsq;
                        //cross[2] = 0;
                        a_n_cross_d2n(0) += cross[0]*weights[k];
                        a_n_cross_d2n(1) += cross[1]*weights[k];
                        //a_n_cross_d2n(2) += cross[2]*weights[k];
                      }
                  },bx0,a_n_cross_d2n,weights,a_h);
  // copy lowest slice to all the other slices.
  spreadSlice(a_n_cross_d2n);
}
template<class T, MemType MEM>
void
cubedSphereDetAFlux1(
                     BoxData<T,1,MEM>& a_detAFlux1,
                     T a_h)
{
  PR_TIMERS("cubedSphereDetAFlux1");
  // HDF5Handler h5;
#if QUADRATURE1D == 5
  Array<T,5> weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
#endif
#if QUADRATURE1D == 3
  Array<T,3> weights={1.0/6.0,2.0/3.0,1.0/6.0};
#endif
#if QUADRATURE1D == 2
  Array<T,2> weights={.5,.5};
#endif
  int q1d = QUADRATURE1D;
  Box bx = a_detAFlux1.box();
  Point low0 = bx.low();
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  a_detAFlux1.setToZero();
  forallInPlace_p([q1d]PROTO_LAMBDA
                  (Point a_pt,
                   Var<T,1,MEM>&  a_detaflux1,
                   Array<T,QUADRATURE1D>& weights,
                   T a_h)
                  {
                    a_detaflux1(0) = 0.0;
                    T alpha = -.25*M_PI + .5*M_PI*a_pt[1]*a_h;
                    T tana = tan(alpha);
                    for (int k = 0; k < q1d;k++)
                      {                      
                        T beta = -.25*M_PI + .5*M_PI*(a_pt[2]*a_h + k*a_h/(q1d-1));                           T tanb = tan(beta);
                        T secb =  1.0/cos(beta);
                        T fluxpt = tana*secb*secb/
                          (sqrt(tana*tana + tanb*tanb + 1.0)*
                           (1.0 + tanb*tanb));        
                        a_detaflux1(0) += fluxpt*weights[k];
                      }                   
                  },bx0,a_detAFlux1,weights,a_h);

  // copy lowest slice to all the other slices.
  spreadSlice(a_detAFlux1);
}
template<class T, MemType MEM>
void
cubedSphereDetAFlux2(
                     BoxData<T,1,MEM>& a_detAFlux2,
                     T a_h)
{
  PR_TIMERS("cubedSphereDetAFlux2");
#if QUADRATURE1D == 5
  Array<T,5> weights={7.0/90.0,32.0/90.0,12.0/90.0,32.0/90.0,7.0/90.0};
#endif
#if QUADRATURE1D == 3
  Array<T,3> weights={1.0/6.0,2.0/3.0,1.0/6.0};
#endif
#if QUADRATURE1D == 2
  Array<T,2> weights={.5,.5};
#endif
  int q1d = QUADRATURE1D;
  Box bx = a_detAFlux2.box();
  Point low0 = bx.low();
  Point high0 = bx.high();
  high0[0] = low0[0];
  Box bx0(low0,high0);
  forallInPlace_p([q1d]PROTO_LAMBDA
                  (Point a_pt,
                   Var<T,1,MEM>&  a_detaflux2,
                   Array<T,QUADRATURE1D>& weights,
                   T a_h)
                  {
                    a_detaflux2(0) = 0.0;
                    T beta = -.25*M_PI + .5*M_PI*a_pt[2]*a_h;
                    T tanb = tan(beta);
                    for (int k = 0; k < q1d;k++)
                      {
                        T alpha = -.25*M_PI + .5*M_PI*(a_pt[1]*a_h + k*a_h/(q1d-1));
                        T tana = tan(alpha);
                        T seca =  1.0/cos(alpha);
                        T fluxpt = tanb*seca*seca/
                          sqrt(tana*tana + tanb*tanb + 1.0)/(1.0 + tana*tana);
                        a_detaflux2(0) += fluxpt*weights[k];                
                      }
                  },bx0,a_detAFlux2,weights,a_h);
  // copy lowest slice to all the other slices.
   spreadSlice(a_detAFlux2);
}
template<class T, MemType MEM>
void
cubedSphereGeometry(BoxData<T, DIM, MEM>& a_X,
                    FluxBoxData<T,DIM,MEM>& a_NT,
                    BoxData<T,1,MEM>& a_J,
                    BoxData<T, 1, MEM>& a_radius,
                    const Box& a_bx,
                    T a_h,
                    int a_block)
{
    //Array<BoxData<T, DIM, MEM>, DIM> NT3;
    FluxBoxData<T,DIM,MEM> NT3;
    BoxData<T, DIM, MEM> X3;
    cubedSphereGeometry(X3,NT3,a_J,a_radius,a_bx,a_h); 
    a_X.define(X3.box());
    // Shuffle data to get correct orientation.
    Array<Array<unsigned int,DIM>, 6>
        permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
    Array<Array<int,DIM> , 6>
        Xsigns = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}};
    //Array<Array<unsigned int,DIM>, 6>
    //   permuteT = {{2,0,1},{2,0,1},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
    //Array<Array<int,DIM> , 6>
    //  XsignsT = {{-1,1,-1},{1,1,1},{-1,1,1},{1,1,1},{1,-1,1},{-1,-1,1}};
    forallInPlace([]PROTO_LAMBDA(
                Var<T,DIM,MEM>& a_XTrans,                      
                Var<T,DIM,MEM>& a_X3,
                Array<unsigned int,DIM> a_permute,
                Array<int,DIM> a_signs)         
    {       
        for (int dir = 0; dir < DIM ; dir++)
        {
            a_XTrans(dir) = a_signs[dir]*a_X3(a_permute[dir]);
        }
    },a_X,X3,permute[a_block],Xsigns[a_block]);
    a_NT.define(NT3.box());
    for (int dir = 0; dir < DIM; dir++)
    {
        forallInPlace([]PROTO_LAMBDA(
                    Var<T,DIM,MEM>& a_ntface,                      
                    Var<T,DIM,MEM>& a_ntface3,
                    Array<unsigned int,DIM> a_permute,
                    Array<int,DIM> a_signs)         
        {       
            for (int dir = 0; dir < DIM; dir++)
            {
                a_ntface(dir) = a_signs[dir]*a_ntface3(a_permute[dir]);
            }
        },a_NT[dir],NT3[dir],permute[a_block],Xsigns[a_block]);
    }
}

// Compute cubed-sphere geometry for the case of radial normal at block center = \f$ \hat{x} \f$ .
template<class T, MemType MEM>
void
cubedSphereGeometry(BoxData<T, DIM, MEM>& a_X,
                    FluxBoxData<T,DIM,MEM>& a_NT,
                    BoxData<T,1,MEM>& a_J,
                    BoxData<T, 1, MEM>& a_radius,
                    const Box& a_bx,
                    T a_h)
{
  PR_TIMERS("cubedSphereGeometry");
  a_X.define(a_bx.extrude(Point::Ones()));
  cubedSphere_nSphere(a_X,a_h);
  forallInPlace([]PROTO_LAMBDA(
                               Var<T,DIM,MEM>& a_x,
                               Var<T,1,MEM>& a_r
                               )
                {
                  a_x(0) *= a_r(0);
                  a_x(1) *= a_r(0);
                  a_x(2) *= a_r(0);
                },a_X,a_radius);

  a_NT.define(a_bx);
  BoxData<T,DIM,MEM> d1n_cross_n(a_bx.extrude(Point::Ones() - Point::Basis(1)));
  BoxData<T,DIM,MEM> n_cross_d2n(a_bx.extrude(Point::Ones() - Point::Basis(2)));
  cubedSphere_edge1(d1n_cross_n,a_h);
  cubedSphere_edge2(n_cross_d2n,a_h);
  // Cofactors.

  Stencil<T> div1 = (1.0)*(Shift(Point::Basis(1)))+ (-1.0)*(Shift(Point::Zeros()));
  Stencil<T> div2 = (1.0)*(Shift(Point::Basis(2)))+ (-1.0)*(Shift(Point::Zeros()));
  
  a_NT[0].setToZero(); 
  a_NT[0] += div2(d1n_cross_n,.5/a_h);
  a_NT[0] += div1(n_cross_d2n,.5/a_h);
  
  forallInPlace([]PROTO_LAMBDA(
              Var<T,DIM,MEM>& a_NT0,
              const Var<T,1,MEM>& a_r)
      {
          T rsq = a_r(0)*a_r(0);
          a_NT0(0) *= rsq;
          a_NT0(1) *= rsq;
          a_NT0(2) *= rsq;  
      },a_NT[0],a_radius);
                
  n_cross_d2n.copyTo(a_NT[1]);
  d1n_cross_n.copyTo(a_NT[2]);
  a_NT[1]*= -1.0/a_h;
  a_NT[2]*= -1.0/a_h;
 
  {
    const BoxData<T,1,MEM> rtop = alias(a_radius,Point::Basis(0,-1));
    forallInPlace([]PROTO_LAMBDA(
                Var<T,DIM,MEM>& a_nttan1,  
                const Var<T,1,MEM>& a_rtop,
                const Var<T,1,MEM>& a_r)
        {
            T drsq = (a_rtop(0)*a_rtop(0) - a_r(0)*a_r(0))/2;
            for (int dir = 0; dir < DIM ; dir++)
            {
                a_nttan1(dir) *= drsq;              
            }
        },a_NT[1],rtop,a_radius);
    forallInPlace([]PROTO_LAMBDA(
            Var<T,DIM,MEM>& a_nttan2,
            const Var<T,1,MEM>& a_rtop,
            const Var<T,1,MEM>& a_r)
    {
        T drsq = (a_rtop(0)*a_rtop(0) - a_r(0)*a_r(0))/2;
        for (int dir = 0; dir < DIM ; dir++)
        {
            a_nttan2(dir) *= drsq;
        }
    },a_NT[2],rtop,a_radius);
  }
   //Jacobian.
  a_J.define(a_bx);
  a_J.setToZero();
  BoxData<T,MEM> detAFlux1(a_bx.extrude(Point::Basis(1)));
  BoxData<T,MEM> detAFlux2(a_bx.extrude(Point::Basis(2)));
  cubedSphereDetAFlux1(detAFlux1,a_h);
  cubedSphereDetAFlux2(detAFlux2,a_h);

  a_J += div1(detAFlux1,.5);
  a_J += div2(detAFlux2,.5);
  const BoxData<T,1,MEM> rtop = alias(a_radius,Point::Basis(0,-1));
  forallInPlace([]PROTO_LAMBDA(
                               Var<T,1,MEM>& a_j,
                               const Var<T,1,MEM>& a_rtop,
                               const Var<T,1,MEM>& a_r)
                {
                  T factor = (a_rtop(0)*a_rtop(0)*a_rtop(0)
                              - a_r(0)*a_r(0)*a_r(0))/3;
                  a_j(0) *= factor;
                },a_J,rtop,a_radius);
  a_J*= 1.0/(a_h*a_h);
}

template<class T, MemType MEM>
BoxData<T,DIM,MEM,DIM>
cofactorMatrix(
        //const Array<BoxData<T,DIM,MEM>, DIM>& a_NTFaces,
        const FluxBoxData<T,DIM,MEM>& a_NTFaces,
        int a_dirNorm,
        Array<T,DIM> a_dx)                         
{
    auto N = cofactorMatrix(a_NTFaces, a_dirNorm);
    T area = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        if (d == a_dirNorm) {continue; }
        area *= a_dx[d];
    }
    N /= area;
    return N;
}

template<class T, MemType MEM>
BoxData<T,DIM,MEM,DIM>
cofactorMatrix(
        //const Array<BoxData<T,DIM,MEM>, DIM>& a_NTFaces,
        const FluxBoxData<T,DIM,MEM>& a_NTFaces,
        int a_dirNorm)                         
{
    PR_TIME("Operator::cofactorMatrix");
    Array<BoxData<T,DIM,MEM>, DIM> NTAveToFace;
    for (int dir = 0; dir < DIM; dir++)
    {
        auto S = Stencil<T>::faceToCell(dir,4)*Stencil<T>::CellToFace(a_dirNorm);
        if (dir != a_dirNorm)
        {
            NTAveToFace[dir] = S(a_NTFaces[dir]);
        }
    }

#if DIM==2
    auto retval = forall<T,DIM,MEM,DIM>
        ([ ] PROTO_LAMBDA(
              Var<T,DIM,MEM,DIM>& a_retval,
              const Var<T,DIM,MEM>& a_ntNorm,
              const Var<T,DIM,MEM>& a_ntTan,
              int a_dirNorm)
         {
             for (int dir2 = 0; dir2 < DIM;dir2++)
             {
                 a_retval(dir2,a_dirNorm) = a_ntNorm(dir2);
                 a_retval(dir2,(a_dirNorm+1)%DIM) = a_ntTan(dir2);
             }
         },
         a_NTFaces[a_dirNorm],
         NTAveToFace[(a_dirNorm+1)%DIM],
         a_dirNorm);
#elif DIM==3
    auto retval =
        forall<T,DIM,MEM,DIM,1>
        ([ ] PROTO_LAMBDA(
              Var<T,DIM,MEM,DIM>& a_retval,
              const Var<T,DIM,MEM>& a_ntNorm,
              const Var<T,DIM,MEM>& a_ntTan1,
              const Var<T,DIM,MEM>& a_ntTan2,
              int a_dirNorm)
         {
             for (int dir1 = 0; dir1 < DIM;dir1++)
             {
                 a_retval(dir1,a_dirNorm) = a_ntNorm(dir1);
                 a_retval(dir1,(a_dirNorm+1)%DIM) = a_ntTan1(dir1);
                 a_retval(dir1,(a_dirNorm+2)%DIM) = a_ntTan2(dir1);
             }
         },
         a_NTFaces[a_dirNorm],
         NTAveToFace[(a_dirNorm+1)%DIM],
         NTAveToFace[(a_dirNorm+2)%DIM],
         a_dirNorm);
#endif
    return retval;
}

//==========================================================================================
// JACOBIAN

template <class T, MemType MEM>
void
jacobian(
        BoxData<T,1,MEM>&         a_J,
        const BoxData<T,DIM,MEM>& a_X,
        const FluxBoxData<T,DIM,MEM>& a_NT)
        //const Array<BoxData<T,DIM,MEM>, DIM>& a_NT)
{
    MayDay<void>::Abort("Operator::jacobian | Error: Not implemented.");
}

template <class T, MemType MEM>
BoxData<T,1,MEM>
jacobian(
        const BoxData<T,DIM,MEM>& a_X,
        const FluxBoxData<T,DIM,MEM>& a_NT,
        Array<T,DIM> a_dx)
{
    T volume = a_dx.dot(Array<T,DIM>::Ones());
    auto J = jacobian(a_X, a_NT);
    J /= volume;
    return J;
}
template <class T, MemType MEM>
BoxData<T,1,MEM>
jacobian(
        const BoxData<T,DIM,MEM>& a_X,
        const FluxBoxData<T,DIM,MEM>& a_NT)
        //const Array<BoxData<T,DIM,MEM>, DIM>& a_NT)
{
    PR_TIME("Operator::jacobian");
    Array<BoxData<T,1,MEM>,DIM> dfdx;
    for (int dir = 0; dir < DIM; dir++)
    {
        BoxData<T,DIM,MEM> xFace4 = Stencil<T>::cornersToFaces(dir,4)(a_X);
        BoxData<T,DIM,MEM> xFace2 = Stencil<T>::cornersToFaces(dir,2)(a_X);
        auto flux = Operator::_faceMatrixProductATB(a_NT[dir],xFace4,a_NT[dir],xFace2,dir);
        dfdx[dir] = Stencil<T>::FluxDivergence(dir)(flux,+1);
    }
#if DIM==3
    auto jac =
        forall<T,1,MEM>([]PROTO_LAMBDA(Var<T,1,MEM>& a_retval,
                    Var<T,1,MEM>& a_f0,
                    Var<T,1,MEM>& a_f1,
                    Var<T,1,MEM>& a_f2)
                {
                a_retval(0) =  (a_f0(0) + a_f1(0) + a_f2(0))/DIM;

                },dfdx[0],dfdx[1],dfdx[2]);
    return jac;
#elif DIM==2
    BoxData<T,1,MEM> jac =
      forall<T,1,MEM,1,1>([]PROTO_LAMBDA(Var<T,1,MEM>& a_retval,
                    Var<T,1,MEM>& a_f0,
                    Var<T,1,MEM>& a_f1)
                {
                a_retval(0) = (1.0/DIM)*(a_f0(0) + a_f1(0));

                },dfdx[0],dfdx[1]);
    return jac;
#else
    std::cout << "jacobian only implemented for DIM = 2,3" << endl;
    abort();
#endif
}

//==========================================================================================
// EDGE CROSS PRODUCT

template<class T, MemType MEM>
BoxData<T,3,MEM>
_edgeCrossProduct3D(
                             const BoxData<T,3,MEM>& a_ql4,
                             const BoxData<T,3,MEM>& a_qr4,
                             const BoxData<T,3,MEM>& a_ql2,
                             const BoxData<T,3,MEM>& a_qr2,
                             int a_edgedir)
{
    PR_TIME("Operator::edgeCrossProduct");
    Array<BoxData<T,1,MEM>,3> qlcomp2,qrcomp2,qlcomp4,qrcomp4;  
    for (int comp = 0; comp < 3; comp++)
    {
        qlcomp2[comp] = slice(a_ql2,comp);
        qrcomp2[comp] = slice(a_qr2,comp);
        qlcomp4[comp] = slice(a_ql4,comp);
        qrcomp4[comp] = slice(a_qr4,comp);
    }
    auto prod_01 = _product1D(qlcomp4[0],qrcomp4[1],qlcomp2[0],qrcomp2[1],a_edgedir);
    auto prod_10 = _product1D(qlcomp4[1],qrcomp4[0],qlcomp2[1],qrcomp2[0],a_edgedir);
    auto prod_12 = _product1D(qlcomp4[1],qrcomp4[2],qlcomp2[1],qrcomp2[2],a_edgedir);
    auto prod_21 = _product1D(qlcomp4[2],qrcomp4[1],qlcomp2[2],qrcomp2[1],a_edgedir);
    auto prod_20 = _product1D(qlcomp4[2],qrcomp4[0],qlcomp2[2],qrcomp2[0],a_edgedir);
    auto prod_02 = _product1D(qlcomp4[0],qrcomp4[2],qlcomp2[0],qrcomp2[2],a_edgedir);

    auto retBD =
        forall<T,3,MEM>([]PROTO_LAMBDA(Var<T,3,MEM>& a_retval,
                    Var<T,1,MEM>& a_pr_01,Var<T,1,MEM>& a_pr_10,
                    Var<T,1,MEM>& a_pr_12,Var<T,1,MEM>& a_pr_21,
                    Var<T,1,MEM>& a_pr_20,Var<T,1,MEM>& a_pr_02)
            {
                a_retval(0) = a_pr_12(0) - a_pr_21(0);
                a_retval(1) = a_pr_20(0) - a_pr_02(0);
                a_retval(2) = a_pr_01(0) - a_pr_10(0);
            },prod_01,prod_10,prod_12,prod_21,prod_20,prod_02);
    return retBD;
}     

template<class T, MemType MEM>
void
_edgeCrossProduct3D(
        BoxData<T,3,MEM>&       a_product,
        const BoxData<T,3,MEM>& a_ql4,
        const BoxData<T,3,MEM>& a_qr4,
        const BoxData<T,3,MEM>& a_ql2,
        const BoxData<T,3,MEM>& a_qr2,
        int a_dir)
{
    MayDay<void>::Abort("Operator::_edgeCrossProduct3D | Error: Not implemented.");
}

template<class T, MemType MEM>
BoxData<T,3,MEM>
edgeCrossProduct3D(
        const BoxData<T,3,MEM>& a_ql,
        const BoxData<T,3,MEM>& a_qr,
        int a_edgedir)
{
    MayDay<void>::Abort("Operator::edgeCrossProduct3D | Error: Not implemented.");
}

template<class T, MemType MEM>
void
edgeCrossProduct3D(
        BoxData<T,3,MEM>&       a_product,
        const BoxData<T,3,MEM>& a_ql,
        const BoxData<T,3,MEM>& a_qr,
        int a_edgedir)
{
    MayDay<void>::Abort("Operator::edgeCrossProduct3D | Error: Not implemented.");
}

//==========================================================================================
// CELL MATRIX PRODUCT
template<class T, unsigned int CL,unsigned int DL,
         unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,CL,MEM,DR,E>
_cellMatrixProductAB(
                     const BoxData<T,CL,MEM,DL,E>& a_ql4,
                     const BoxData<T,CR,MEM,DR,E>& a_qr4,
                     const BoxData<T,CL,MEM,DL,E>& a_ql2,
                     const BoxData<T,CR,MEM,DR,E>& a_qr2)
{
    PR_TIME("Operator::matrixProductAB");
    PR_assert(DL == CR);
    Box bx = a_ql4.box()&a_qr4.box();
    for (int dir = 0; dir < DIM; dir++)
    {
        bx = bx&a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
    }
    BoxData<T,CL,MEM,DR,E> prod(bx);
    prod.setToZero();
    for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
    for (int j = 0; j < DR; j++)
    for (int k = 0; k < DL;k++)
    {
        auto prodcomp = slice(prod,i,j,ee);
        auto qr4comp = slice(a_qr4,k,j,ee);
        auto qr2comp = slice(a_qr2,k,j,ee);
        auto ql4comp = slice(a_ql4,i,k,ee);
        auto ql2comp = slice(a_ql2,i,k,ee);
        prodcomp+=Operator::_cellProduct(ql4comp,qr4comp,ql2comp,qr2comp);
    }
    return prod;
}

//==========================================================================================
// FACE MATRIX PRODUCT
template<class T, unsigned int CL,unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,CL,MEM,DR,E>
_faceMatrixProductAB(
        const BoxData<T,CL,MEM,DL,E>& a_ql4,
        const BoxData<T,CR,MEM,DR,E>& a_qr4,
        const BoxData<T,CL,MEM,DL,E>& a_ql2,
        const BoxData<T,CR,MEM,DR,E>& a_qr2,
        int a_dir)
{
    PR_TIME("Operator::faceMatrixProductAB");
    PR_assert(uint(DL) == CR);
    Box bx = a_ql4.box() & a_qr4.box();
    for (int dir = 0; dir < DIM; dir++)
    {
        if (dir!=a_dir)
        {
            bx &= a_ql2.box().grow(dir,-1);
            bx &= a_qr2.box().grow(dir,-1);
        }
    }
    BoxData<T,CL,MEM,DR,E> prod(bx);
    prod.setToZero();
    for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
    for (int j = 0; j < DR; j++)
    for (int k = 0; k < DL;k++)
    {
        auto prodcomp = slice(prod,i,j,ee);
        auto qr4comp = slice(a_qr4,k,j,ee);
        auto qr2comp = slice(a_qr2,k,j,ee);
        auto ql4comp = slice(a_ql4,i,k,ee);
        auto ql2comp = slice(a_ql2,i,k,ee);
        prodcomp+=Operator::_faceProduct(ql4comp,qr4comp,ql2comp,qr2comp,a_dir);
    }
    return prod;
}
template<class T, unsigned int CL,unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
void 
matrixProductATB(
                 BoxData<T,DL,MEM,DR,E>& a_qprod,
                 const BoxData<T,DL,MEM,CL,E>& a_ql4,
                 const BoxData<T,CR,MEM,DR,E>& a_qr4)
{
  PR_TIME("Operator::matrixProductATB");
  PR_assert(CL == CR);
  a_qprod.setToZero();
  for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
      for (int j = 0; j < DR; j++)
        for (int k = 0; k < DL;k++)
          {
            auto prodcomp = slice(a_qprod,i,j,ee);
            auto qr4comp = slice(a_qr4,k,j,ee);              
            auto ql4comp = slice(a_ql4,k,i,ee);              
            prodcomp+=Operator::cellProduct(ql4comp,qr4comp);
          }
}
template<class T, unsigned int CL,unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
void 
matrixProductAB(
                BoxData<T,CL,MEM,DR,E>& a_qprod,
                const BoxData<T,DL,MEM,CL,E>& a_ql4,
                const BoxData<T,CR,MEM,DR,E>& a_qr4)
{
  PR_TIME("Operator::matrixProductATB");
  PR_assert(CL == CR);
  a_qprod.setToZero();
  for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
      for (int j = 0; j < DR; j++)
        for (int k = 0; k < DL;k++)
          {
            auto prodcomp = slice(a_qprod,i,j,ee);
            auto qr4comp = slice(a_qr4,k,j,ee);              
            auto ql4comp = slice(a_ql4,i,k,ee);              
            prodcomp+=Operator::cellProduct(ql4comp,qr4comp);
          }
}
template<class T, unsigned int CL,unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
void
_faceMatrixProductAB(
        BoxData<T,CL,MEM,DR,E>&       a_product,
        const BoxData<T,CL,MEM,DL,E>& a_ql4,
        const BoxData<T,CR,MEM,DR,E>& a_qr4,
        const BoxData<T,CL,MEM,DL,E>& a_ql2,
        const BoxData<T,CR,MEM,DR,E>& a_qr2,
        int a_dir)
{
    MayDay<void>::Abort("Operator::_faceMatrixProduct | Error: Not implemented.");
}

template<class T, unsigned int CL,unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,CL,MEM,DR,E>
faceMatrixProductAB(
        const BoxData<T,CL,MEM,DL,E>& a_ql,
        const BoxData<T,CR,MEM,DR,E>& a_qr,
        int a_dir)
{
    MayDay<void>::Abort("Operator::faceMatrixProduct | Error: Not implemented.");
}

template<class T, unsigned int CL,unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
void
faceMatrixProductAB(
        BoxData<T,CL,MEM,DR,E>&       a_product,
        const BoxData<T,CL,MEM,DL,E>& a_ql,
        const BoxData<T,CR,MEM,DR,E>& a_qr,
        int a_dir)
{
    MayDay<void>::Abort("Operator::faceMatrixProduct | Error: Not implemented.");
}

template<class T, unsigned int CL, unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,DL,MEM,DR,E>
_faceMatrixProductATB(
        const BoxData<T,CL,MEM,DL,E>& a_ql4,
        const BoxData<T,CR,MEM,DR,E>& a_qr4,
        const BoxData<T,CL,MEM,DL,E>& a_ql2,
        const BoxData<T,CR,MEM,DR,E>& a_qr2,
        int a_dir)
{
    PR_TIME("Operator::faceMatrixProductATB");
    PR_assert(CR == CL);
    Box bx = a_ql4.box()&a_qr4.box();
    for (int dir = 0; dir < DIM; dir++)
    {
        if (dir!=a_dir) bx = bx
            &a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
    }
    BoxData<T,DL,MEM,DR,E> prod(bx);
    prod.setToZero();
    for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < DL; i++)
    for (int j = 0; j < DR; j++)
    for (int k = 0; k < CL;k++)
    {
        auto prodcomp = slice(prod,i,j,ee);
        auto qr4comp = slice(a_qr4,k,j,ee);
        auto qr2comp = slice(a_qr2,k,j,ee);
        auto ql4comp = slice(a_ql4,k,i,ee);
        auto ql2comp = slice(a_ql2,k,i,ee);
        prodcomp+=Operator::_faceProduct(ql4comp,qr4comp,ql2comp,qr2comp,a_dir);
    }
    return prod;
}

template<class T, unsigned int CL, unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
void
_faceMatrixProductATB(
        BoxData<T,DL,MEM,DR,E>&       a_product,
        const BoxData<T,CL,MEM,DL,E>& a_ql4,
        const BoxData<T,CR,MEM,DR,E>& a_qr4,
        const BoxData<T,CL,MEM,DL,E>& a_ql2,
        const BoxData<T,CR,MEM,DR,E>& a_qr2,
        int a_dir)
{
    PR_TIME("Operator::faceMatrixProductATB");
    PR_assert(CR == CL);
    Box bx = a_ql4.box()&a_qr4.box();
    for (int dir = 0; dir < DIM; dir++)
    {
        if (dir!=a_dir) bx = bx
            &a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
    }
    BoxData<T,DL,MEM,DR,E> prod(bx);
    prod.setToZero();
    for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < DL; i++)
    for (int j = 0; j < DR; j++)
    for (int k = 0; k < CL;k++)
    {
        auto prodcomp = slice(prod,i,j,ee);
        auto qr4comp = slice(a_qr4,k,j,ee);
        auto qr2comp = slice(a_qr2,k,j,ee);
        auto ql4comp = slice(a_ql4,k,i,ee);
        auto ql2comp = slice(a_ql2,k,i,ee);
        prodcomp+=Operator::_faceProduct(ql4comp,qr4comp,ql2comp,qr2comp,a_dir);
    }
}

template<class T, unsigned int CL, unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,CL,MEM,CR,E>
_faceMatrixProductABT(
        const BoxData<T,CL,MEM,DL,E>& a_ql4,
        const BoxData<T,CR,MEM,DR,E>& a_qr4,
        const BoxData<T,CL,MEM,DL,E>& a_ql2,
        const BoxData<T,CR,MEM,DR,E>& a_qr2,
        int a_dir)
{
    PR_TIME("Operator::faceMatrixProductABT");
    PR_assert(DL == DR);
    Box bx = a_ql4.box()&a_qr4.box();
    for (int dir = 0; dir < DIM; dir++)
    {
        if (dir!=a_dir) bx = bx
            &a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
    }
    BoxData<T,CL,MEM,CR,E> prod(bx);
    prod.setToZero();
    for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
    for (int j = 0; j < CR; j++)
    for (int k = 0; k < DL;k++)
    {
        auto prodcomp = slice(prod,i,j,ee);
        auto qr4comp = slice(a_qr4,j,k,ee);
        auto qr2comp = slice(a_qr2,j,k,ee);
        auto ql4comp = slice(a_ql4,i,k,ee);
        auto ql2comp = slice(a_ql2,i,k,ee);
        prodcomp+=Operator::_faceProduct(ql4comp,qr4comp,ql2comp,qr2comp,a_dir);
    }
    return prod;
}

template<class T, unsigned int CL, unsigned int DL,
    unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
void
_faceMatrixProductABT(
        BoxData<T,CL,MEM,CR,E>&       a_product,
        const BoxData<T,CL,MEM,DL,E>& a_ql4,
        const BoxData<T,CR,MEM,DR,E>& a_qr4,
        const BoxData<T,CL,MEM,DL,E>& a_ql2,
        const BoxData<T,CR,MEM,DR,E>& a_qr2,
        int a_dir)
{
    MayDay<void>::Abort("_faceMatrixProductABT | Error: Not Implemented.");
}

//==========================================================================================
// PRODUCT 1D
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BoxData<T,C,MEM,D,E>
_product1D(
        const BoxData<T,C,MEM,D,E>& a_ql4,
        const BoxData<T,C,MEM,D,E>& a_qr4,
        const BoxData<T,C,MEM,D,E>& a_ql2,
        const BoxData<T,C,MEM,D,E>& a_qr2,
        int a_dir)
{
    PR_TIME("Operator::product1D");
    EdgeAverageProductOp<T> op(a_dir);
    return op(a_ql4, a_qr4, a_ql2, a_qr2);
}

template<class T, unsigned int C, unsigned int D, unsigned int E, MemType MEM>
BoxData<T,C,MEM,D,E>
_cellTensorQuotient(
                              const BoxData<T,C,MEM,D,E>& a_ql4,
                              const BoxData<T,1,MEM,1,1>& a_qr4,
                              const BoxData<T,C,MEM,D,E>& a_ql2,
                              const BoxData<T,1,MEM,1,1>& a_qr2)
{
  Box bx = a_ql4.box()&a_qr4.box()&a_ql2.box().grow(-1)&a_qr2.box().grow(-1);
  BoxData<T,C,MEM,D,E> quot(bx);
  for (int ee = 0; ee < E;ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
        {
          auto quotcomp = slice(quot,cc,dd,ee);
          auto ql4comp = slice(a_ql4,cc,dd,ee);
          auto ql2comp = slice(a_ql2,cc,dd,ee);
          Operator::_cellQuotient(quotcomp,ql4comp,a_qr4,ql2comp,a_qr2);
        }
  return quot;
}
template<class T, unsigned int CL,unsigned int DL,
           unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,CL,MEM,DR,E>
_matrixProductAB2(
                  const BoxData<T,CL,MEM,DL,E>& a_ql,
                  const BoxData<T,CR,MEM,DR,E>& a_qr)
{
  PR_TIME("Operator::matrixProductAB2");
  PR_assert(uint(DL) == CR);
  Box bx = a_ql.box()&a_qr.box();
  BoxData<T,CL,MEM,DR,E> prod(bx);
  prod.setToZero();
  for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
      for (int j = 0; j < DR; j++)
        for (int k = 0; k < DL;k++)
        {
          auto prodcomp = slice(prod,i,j,ee);
          auto qrcomp = slice(a_qr,k,j,ee);
          auto qlcomp = slice(a_ql,i,k,ee);
          prodcomp+=qrcomp*qlcomp;
        }
  return prod;
}
template<class T, unsigned int CL, unsigned int DL,
         unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,DL,MEM,DR,E>
_matrixProductATB2(
                   const BoxData<T,CL,MEM,DL,E>& a_ql,
                   const BoxData<T,CR,MEM,DR,E>& a_qr)
{
  PR_TIME("Operator::matrixProductATB2");
  PR_assert(CR == CL);
  Box bx = a_ql.box()&a_qr.box();
  BoxData<T,DL,MEM,DR,E> prod(bx);
  prod.setToZero();
  for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < DL; i++)
      for (int j = 0; j < DR; j++)
        for (int k = 0; k < CL;k++)
          {
            auto prodcomp = slice(prod,i,j,ee);
            auto qrcomp = slice(a_qr,k,j,ee);
            auto qlcomp = slice(a_ql,k,i,ee);
            prodcomp+=qrcomp*qlcomp;
          }
  return prod;
}
template<class T, unsigned int CL, unsigned int DL,
         unsigned int CR, unsigned int DR, unsigned int E, MemType MEM>
BoxData<T,CL,MEM,CR,E>
_matrixProductABT2(
                   const BoxData<T,CL,MEM,DL,E>& a_ql,
                   const BoxData<T,CR,MEM,DR,E>& a_qr)
{
  PR_TIME("Operator::faceMatrixProductABT");
  PR_assert(DL == DR);
  Box bx = a_ql.box()&a_qr.box();
  BoxData<T,CL,MEM,CR,E> prod(bx);
  prod.setToZero();
  for (int ee = 0; ee < E;ee++)
    for (int i = 0; i < CL; i++)
      for (int j = 0; j < CR; j++)
        for (int k = 0; k < DL;k++)
          {
            auto prodcomp = slice(prod,i,j,ee);
            auto qrcomp = slice(a_qr,j,k,ee);
            auto qlcomp = slice(a_ql,i,k,ee);            
            prodcomp+=qrcomp*qlcomp;
          }
  return prod;  
}

template<class T,unsigned int C, MemType MEM>
BoxData<T,C,MEM,DIM>
_faceGradient(
        const BoxData<T,C,MEM>& a_phi,
        int a_dirNorm,
        int a_order)
{
    PR_TIME("Operator::faceGradient");

    FaceAverageDiffOp<T> diffOp(a_dirNorm, a_order);
    return diffOp(a_phi);
}

template<class T,unsigned int CFLUX,unsigned int CPRIM, MemType MEM>
BoxData<T,CFLUX,MEM>
MHDSphericalFlux(                     
                 const BoxData<T,CPRIM,MEM>& a_prim4,
                 const BoxData<T,CPRIM,MEM>& a_prim2,
                 const BoxData<T,DIM,MEM,DIM>& a_A4,
                 const BoxData<T,DIM,MEM,DIM>& a_A2,
                 const BoxData<T,1,MEM>& a_DrDetA4,                    
                 const BoxData<T,1,MEM>& a_DrDetA2,
                 const BoxData<T,DIM,MEM>& a_DrAdjA4,                    
                 const BoxData<T,DIM,MEM>& a_DrAdjA2,
                 const T& a_gamma,
                 int a_dir)
{
#define CRHO 0
#define CVELSTART 1
#define CBSTART 4
#define CPRES 7
#define CENG 7
  
  auto wnorm4 = slice(a_prim4,CVELSTART+a_dir);
  auto bnorm4 = slice(a_prim4,CBSTART+a_dir);
  auto wnorm2 = slice(a_prim2,CVELSTART+a_dir);
  auto bnorm2 = slice(a_prim2,CBSTART+a_dir);
  // Volumetric flow rates V_u,V_b.
  auto Vu4 = _faceProduct(a_DrDetA4,wnorm4,a_DrDetA2,wnorm2,a_dir);
  auto Vb4 = _faceProduct(a_DrDetA4,bnorm4,a_DrDetA2,bnorm2,a_dir);
  auto Vu2 = a_DrDetA2*wnorm2; // Placeholder for forall.
  auto Vb2 = a_DrDetA2*bnorm2; // Placeholder for forall.

  // Advective fluxes of density, energy.
  
  auto rho4 = slice(a_prim4,CRHO);
  auto rho2 = slice(a_prim2,CRHO);
  auto fluxRho4 = _faceProduct(Vu4,rho4,Vu2,rho2,a_dir);
  auto fluxRho2 = _matrixProductAB2(Vu2,rho2);
  auto p4 = slice(a_prim4,CPRES); 
  auto p2 = slice(a_prim2,CPRES);

  // Cartesian B, velocities.
  auto w4 = slice<T,CPRIM,DIM,MEM>(a_prim4,CVELSTART);
  auto b4 = slice<T,CPRIM,DIM,MEM>(a_prim4,CBSTART);  
  auto w2 = slice<T,CPRIM,DIM,MEM>(a_prim2,CVELSTART);
  auto b2 = slice<T,CPRIM,DIM,MEM>(a_prim2,CBSTART);
  
  BoxData<T,DIM,MEM> Uface4 = _faceMatrixProductAB(a_A4,w4,a_A2,w2,a_dir);
  BoxData<T,DIM,MEM> Bface4 = _faceMatrixProductAB(a_A4,b4,a_A2,b2,a_dir);
  BoxData<T,DIM,MEM> Uface2 = _matrixProductAB2(a_A2,w2);
  BoxData<T,DIM,MEM> Bface2 = _matrixProductAB2(a_A2,b2);

  // Fluxes for velocity, magnetic fields.
  BoxData<T,DIM,MEM> fluxuu = _faceMatrixProductAB(Uface4,fluxRho4,Uface2,fluxRho2,a_dir);
  auto fluxub = _faceMatrixProductAB(Bface4,Vu4,Bface2,Vu2,a_dir);
  auto fluxbu = _faceMatrixProductAB(Bface4,Vb4,Uface2,Vb2,a_dir);
  auto fluxbb = _faceMatrixProductAB(Bface4,Vb4,Bface2,Vb2,a_dir);
  fluxbb *= 1.0/(4*M_PI);
  
  // Energy as a function of the primitive (Cartesian) variables.
  auto Beng4 = _faceMatrixProductATB(Bface4,Bface4,Bface2,Bface2,a_dir);
  auto Beng2 = _matrixProductATB2(Bface2,Bface2);  
  Beng4 *= 1.0/(8.0*M_PI);
  Beng2 *= 1.0/(8.0*M_PI);
  
  auto Ueng4 = _faceMatrixProductATB(Uface4,Uface4,Uface2,Uface2,a_dir);
  auto Ueng2 = _matrixProductATB2(Uface2,Uface2);
  Ueng4 *= .5;
  Ueng2 *= .5; 
  
  // Kinetic energy advective contribution.
  auto fluxUEng = _faceProduct(fluxRho4,Ueng4,fluxRho2,Ueng2,a_dir);
  
  // pzero = sum of thermal and magnetic pressures.
  auto pzero4 = p4 + Beng4;
  auto pzero2 = p2 + Beng2;

  // Advective + p d(1/rho) work contribution to the energy flux.
  p4 *= 1./(a_gamma - 1.0);
  p2 *= 1./(a_gamma - 1.0);
  auto thermBEng4 = Beng4 + p4 + pzero4;
  auto thermBEng2 = Beng2 + p2 + pzero2;
  
  auto fluxThermBAdv = _faceProduct(Vu4,thermBEng4,Vu2,thermBEng2,a_dir);
  
  // Non-gradient magnetic field contribution to the energy flux.
  BoxData<double,1,MEM> UDotB4 = _faceMatrixProductATB(Uface4,Bface4,Uface2,Bface2,a_dir);
  BoxData<double,1,MEM> UDotB2 = _matrixProductATB2(Uface2,Bface2);
  auto fluxBEng = _faceProduct(Vb4,UDotB4,Vb2,UDotB2,a_dir);
  fluxBEng *= (-1.0/(4.0*M_PI));

  // Pressure forces on the fluid.
  auto pForce = _faceMatrixProductAB(a_DrAdjA4,pzero4,a_DrAdjA2,pzero2,a_dir);
  
  // Assemble into flux vector.
  auto retval = forall<T,CFLUX,MEM,1>
    ([ ] PROTO_LAMBDA(
                      Var<T,CFLUX,MEM,1>& a_retval,
                      Var<T,1,MEM>& a_fluxRho,
                      Var<T,DIM,MEM>& a_fluxuu,
                      Var<T,DIM,MEM>& a_fluxub,
                      Var<T,DIM,MEM>& a_fluxbu,
                      Var<T,DIM,MEM>& a_fluxbb,
                      Var<T,1,MEM>& a_fluxThermBAdv,
                      Var<T,1,MEM>& a_fluxUEng,
                      Var<T,1,MEM>& a_fluxBEng,
                      Var<T,DIM,MEM>& a_pforce)
     {
       a_retval(0) = a_fluxRho(0);
       for (int dir = 0; dir < DIM; dir++)
         {
           a_retval(CVELSTART+dir) = a_fluxuu(dir) + a_pforce(dir) + a_fluxbb(dir);
           a_retval(CBSTART+dir) = a_fluxub(dir) - a_fluxbu(dir);
         }
       a_retval(CENG) = a_fluxThermBAdv(0) + a_fluxUEng(0) + a_fluxBEng(0);
     },
     fluxRho4,fluxuu,fluxub,fluxbu,fluxbb,fluxThermBAdv,fluxUEng,fluxBEng,pForce);
  
  return retval;
}
template<class T,unsigned int C, MemType MEM>
BoxData<T,C,MEM,DIM>
_faceGradxPhi(                     
              const BoxData<T,C,MEM>& a_phi4,
              const BoxData<T,C,MEM>& a_phi2,
              const BoxData<T,DIM,MEM,DIM>& a_NFace4,
              const BoxData<T,DIM,MEM,DIM>& a_NFace2,
             const BoxData<T,1,MEM>& a_JFace4,
              const BoxData<T,1,MEM>& a_JFace2,
              int a_dir)
{
  auto gradxiphi4 = Operator::_faceGradient(a_phi4,a_dir);
  auto gradxiphi2 = Operator::_faceGradient(a_phi2,a_dir,2);
  auto jgradphi4 = _faceMatrixProductABT(gradxiphi4,a_NFace4,gradxiphi2,a_NFace2,a_dir);
  auto jgradphi2 = _matrixProductABT2(gradxiphi2,a_NFace2);
  auto retval = _faceTensorQuotient(jgradphi4,a_JFace4,jgradphi2,a_JFace2,a_dir);
  return retval;
}
