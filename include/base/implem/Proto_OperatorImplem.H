using namespace Proto;

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_avg |= Stencil<T>::Laplacian()(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_convolve(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,MEM,D,E> a_avg(a_ctr.box() & a_2nd.box().grow(-1));
    _convolve(a_avg, a_ctr, a_2nd);
    return a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::convolve(const BoxData<T,C,MEM,D,E>& a_ctr)
{
    return _convolve(a_ctr, a_ctr);
}
/*
template<class T, unsigned int C, MemType MEM, Centering CTR>
void
Operator::_convolve(
        LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_2nd)
{
    //TODO: switch functionality based on value of CTR
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
Operator::convolve(
        LevelBoxData<T,C,MEM,CTR>& a_avg,
        const LevelBoxData<T,C,MEM,CTR>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}
*/
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_ctr |= Stencil<T>::Laplacian()(a_2nd, -1.0/24);
    a_ctr += a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_deconvolve(
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,MEM,D,E> a_ctr(a_avg.box() & a_2nd.box().grow(-1));
    _deconvolve(a_ctr, a_avg, a_2nd);
    return a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::deconvolve(const BoxData<T,C,MEM,D,E>& a_avg)
{
    return _deconvolve(a_avg, a_avg);
}
/*
template<class T, unsigned int C, MemType MEM, Centering CTR>
void
Operator::_deconvolve(
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
Operator::deconvolve(
        LevelBoxData<T,C,MEM,CTR>& a_ctr,
        const LevelBoxData<T,C,MEM,CTR>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}
*/
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_avg |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    _convolveFace(a_avg, a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_convolveFace(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::convolveFace(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    return _convolveFace(a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_ctr |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, -1.0/24.0);
    a_ctr += a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    _deconvolveFace(a_ctr, a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_deconvolveFace(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::deconvolveFace(
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    return _deconvolveFace(a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_cellProduct(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::cellProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellProduct(a_pavg, a_qavg, a_pavg, a_qavg);
}
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellProduct | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellProduct | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellProduct | incompatible patch sizes.");
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);

        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, int v_dir)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_dir == 0)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    }
                }
            }, a_product, Dp, Dq, ii);
    }

    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += v_p(cc,dd,ee)*v_q(cc,dd,ee);
            }
        }, a_product, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellProduct(a_product, a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_faceProduct(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::faceProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceProduct(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    bool do_replace = true;
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){continue;}
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, bool v_replace)
            {
                for (int cc = 0; cc < C; cc++)
                {
                for (int dd = 0; dd < D; dd++)
                {
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_replace)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                    }
                }
                }
                }
            }, a_dest, Dp, Dq, do_replace);
        do_replace = false;
    }
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += v_p(cc,dd,ee)*v_q(cc,dd,ee);
            }
        }, a_dest, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceProduct(a_dest, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}


template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_cellQuotient(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::cellQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellQuotient(a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
  //PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    //BoxData<T,C,MEM,D,E> out(a_p2nd.box().grow(-1));
    //out.setVal(0);
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q, int v_dir)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_dir == 0)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    }
                    v_out(cc,dd,ee) -= v_p(cc,dd,ee)/v_q(cc,dd,ee)*v_dq(cc,dd,ee)*v_dq(cc,dd,ee);
                    v_out(cc,dd,ee) *= -1.0/(12.0*v_q(cc,dd,ee)*v_q(cc,dd,ee));
                }
            }, a_quotient, Dp, Dq, a_p2nd, a_q2nd, ii);
    }
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_faceQuotient(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::faceQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceQuotient(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::_faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){continue;}
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q, int v_dir)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {
                    if (v_dir == 0)
                    {
                        v_out(cc,dd,ee) = v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    } else {
                        v_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee);
                    }
                    v_out(cc,dd,ee) -= v_p(cc,dd,ee)/v_q(cc,dd,ee)*v_dq(cc,dd,ee)*v_dq(cc,dd,ee);
                    v_out(cc,dd,ee) *= -1.0/(12.0*v_q(cc,dd,ee)*v_q(cc,dd,ee));
                }
            }, a_quotient, Dp, Dq, a_p2nd, a_q2nd, ii);
    }
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
Operator::faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}
template <class T, MemType MEM>
  BoxData<T,DIM,MEM>
 Operator::_cofactor(
                    const BoxData<T,DIM,MEM>& a_X,
                    int a_dir)
 {
#if DIM==2
   int sign=1;
   if (a_dir == 1) sign=-1;
   int dirperp = (a_dir+1)%DIM;
   Box bxret = a_X.box()&a_X.box().shift(1,dirperp);
   BoxData<T,DIM,MEM> ret(bxret);
   for (int comp = 0;comp < DIM;comp++)
     {
       int compperp = (comp+1)%DIM;
       auto retslice = slice(ret,compperp);
       auto Xslice = slice(a_X,comp);
       T compsign = -2*comp+1; 
       retslice = Stencil<T>::FluxDivergence(dirperp)(Xslice,compsign*sign);
     }
   return ret;
#elif DIM==3
   int dirperp1 = (a_dir+1)%DIM;
   int dirperp2 = (a_dir+2)%DIM;
   BoxData<T,DIM> XAv1 = Stencil<T>::faceToCell(dirperp1)(a_X);
   BoxData<T,DIM> dXdxi1 = Stencil<T>::divergence(dirperp1)(a_X);
   auto XBydXdxi1 = Operator::edgeCrossProduct3D(XAv1,dXdxi1,XAv1,dXdxi1,dirperp1);
   BoxData<T,DIM> XAv2 = Stencil<T>::faceToCell(dirperp2)(a_X);
   BoxData<T,DIM> dXdxi2 = Stencil<T>::divergence(dirperp2)(a_X);
   auto dXdxi2ByX = Operator::edgeCrossProduct3D(dXdxi2,XAv2,dXdxi2,XAv2,dirperp2);
   auto d12 = Stencil<T>::divergence(dirperp2)(XBydXdxi1,.5);
   auto d21 = Stencil<T>::divergence(dirperp1)(dXdxi2ByX,.5);
   auto ret = forall([]PROTO_LAMBDA(Var<T,DIM>& a_d21,Var<T,DIM>& a_d12)
                     {
                       Var<T,DIM> retvalpt;
                       for (int comp = 0; comp < DIM; comp++)
                         {
                           retval(comp) = d12(comp) + d21(comp);
                         }
                       return retvalpt;
                     },d21,d12);
   return ret;
#else
   std::cout << "cofactor only implemented for DIM = 2,3" << endl;
   abort();
#endif
 }
template <class T, MemType MEM>
BoxData<T,1,MEM>
Operator::_jacobian(
                    const BoxData<T,DIM,MEM>& a_X,
                    const std::array<BoxData<T,DIM,MEM>, DIM>& a_NT)
{
  std::array<BoxData<T,1,MEM>,DIM> dfdx;
  for (int dir = 0; dir < DIM; dir++)
    {
      BoxData<T,DIM> xFace4 = Stencil<T>::cornersToFaces(dir,4)(a_X);
      BoxData<T,DIM> xFace2 = Stencil<T>::cornersToFaces(dir,2)(a_X);
      auto flux = Operator::_faceMatrixProductATB(a_NT[dir],xFace4,a_NT[dir],xFace2,dir);
      dfdx[dir] = Stencil<T>::FluxDivergence(dir)(flux);
    }
#if DIM==3
  auto jac =
    forall<T,1>([]PROTO_LAMBDA(Var<T,1>& a_retval,
                          Var<T,1>& a_f0,
                          Var<T,1>& a_f1,
                          Var<T,1>& a_f2)
           {
             a_retval(0) =  a_f0(0) + a_f1(0) + a_f2(0);
             
           },dfdx[0],dfdx[1],dfdx[2]);
  return jac;
#elif DIM==2
auto jac =
  forall<T,1>([]PROTO_LAMBDA(Var<T,1>& a_retval,
                          Var<T,1>& a_f0,
                          Var<T,1>& a_f1)
           {
             a_retval(0) = a_f0(0) + a_f1(0);
             
           },dfdx[0],dfdx[1]);
  return jac;
#else
  std::cout << "jacobian only implemented for DIM = 2,3" << endl;
  abort();
#endif
}
template<class T, MemType MEM>
BoxData<T,3,MEM>
Operator::_edgeCrossProduct3D(
                             const BoxData<T,3,MEM>& a_ql4,
                             const BoxData<T,3,MEM>& a_qr4,
                             const BoxData<T,3,MEM>& a_ql2,
                             const BoxData<T,3,MEM>& a_qr2,
                             int a_edgedir)
{
  std::array<BoxData<T,1,MEM>,DIM> qlcomp2,qrcomp2,qlcomp4,qrcomp4;  
  for (int comp = 0; comp < DIM; comp++)
    {
      qlcomp2[comp] = slice(a_ql2,comp);
      qrcomp2[comp] = slice(a_qr2,comp);
      qlcomp4[comp] = slice(a_ql4,comp);
      qrcomp4[comp] = slice(a_qr4,comp);
    }
  auto prod_01 = product1D(qlcomp4[0],qrcomp4[1],qlcomp2[0],qrcomp2[1],a_edgedir);
  auto prod_10 = product1D(qlcomp4[1],qrcomp4[0],qlcomp2[1],qrcomp2[0],a_edgedir);
  auto prod_12 = product1D(qlcomp4[1],qrcomp4[2],qlcomp2[1],qrcomp2[2],a_edgedir);
  auto prod_21 = product1D(qlcomp4[2],qrcomp4[1],qlcomp2[2],qrcomp2[1],a_edgedir);
  auto prod_20 = product1D(qlcomp4[2],qrcomp4[0],qlcomp2[2],qrcomp2[0],a_edgedir);
  auto prod_02 = product1D(qlcomp4[0],qrcomp4[2],qlcomp2[0],qrcomp2[2],a_edgedir);

  auto retBD =
    forall([]PROTO_LAMBDA(Var<T,DIM>& a_retval,
                          Var<T,1>& a_pr_01,Var<T,1>& a_pr_10,
                          Var<T,1>& a_pr_12,Var<T,1>& a_pr_21,
                          Var<T,1>& a_pr_20,Var<T,1>& a_pr_02)
           {
             a_retval(0) = a_pr_12(0) - a_pr_21(0);
             a_retval(1) =-(a_pr_20(0) - a_pr_02(0));
             a_retval(2) = a_pr_01(0) - a_pr_10(0);
           },prod_01,prod_10,prod_12,prod_21,prod_20,prod_02);
  return retBD;
}      
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
Operator::_Product1D(
                    const BoxData<T,C,MEM,D,E>& a_ql4,
                    const BoxData<T,C,MEM,D,E>& a_qr4,
                    const BoxData<T,C,MEM,D,E>& a_ql2,
                    const BoxData<T,C,MEM,D,E>& a_qr2,
                    int a_dir)
{
  auto b = (a_ql2.box().grow(a_dir,-1) & a_qr2.box().grow(a_dir,-1));
    b &= (a_ql4.box() & a_qr4.box());
    BoxData<T,C,MEM,D,E> prod(b);
    prod.setToZero();
    auto Diff = 0.5*Shift::Basis(a_dir,1) - 0.5*Shift::Basis(a_dir,-1);
    BoxData<T,C,MEM,D,E> Dql = Diff(a_ql2);
    BoxData<T,C,MEM,D,E> Dqr = Diff(a_qr2);
    
    forallInPlace(
                  [ ] PROTO_LAMBDA (Var<T,C,MEM,D,E>& prod_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq)
                  {
                    for (int cc = 0; cc < C; cc++)
                      for (int dd = 0; dd < D; dd++)
                        for (int ee = 0; ee < E; ee++)
                          {
                            prod_out(cc,dd,ee) += v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/12.0;
                          }
                  }, prod, Dql, Dqr);
    forallInPlace(
                  [ ] PROTO_LAMBDA (Var<T,C,MEM,D,E>& prod_out, Var<T,C,MEM,D,E>& v_l, Var<T,C,MEM,D,E>& v_r)
                  {
                    for (int cc = 0; cc < C; cc++)
                      for (int dd = 0; dd < D; dd++)
                        for (int ee = 0; ee < E; ee++)
                          {
                            prod_out(cc,dd,ee) += v_l(cc,dd,ee)*v_r(cc,dd,ee);
                          }
                  }, prod, a_ql4, a_qr4);
    return prod;
}
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E>
Operator::_cellTensorQuotient(
                              const BoxData<T,1,MEM,1,1>& a_ql4,
                         const BoxData<T,C,MEM,D,E>& a_qr4,
                         const BoxData<T,1,MEM,1,1>& a_ql2,
                         const BoxData<T,C,MEM,D,E>& a_qr2)
{
  Box bx = a_ql4.box()&a_qr4.box()&a_ql2.box().grow(-1)&a_qr2.box().grow(-1);
  BoxData<T,C,MEM,D,E> quot(bx);
  for (int ee = 0; ee < E;ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
        {
          auto quotcomp = slice(quot,cc,dd,ee);
          auto qr4comp = slice(a_qr4,cc,dd,ee);
          auto qr2comp = slice(a_qr2,cc,dd,ee);
          Operator::_cellQuotient(quotcomp,a_ql4,qr4comp,a_ql2,qr2comp);
        }
  return quot;
}
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E>
Operator::_faceTensorQuotient(
                              const BoxData<T,1,MEM,1,1>& a_ql4,
                              const BoxData<T,C,MEM,D,E>& a_qr4,
                              const BoxData<T,1,MEM,1,1>& a_ql2,
                              const BoxData<T,C,MEM,D,E>& a_qr2,
                              int a_dir)
{
  Box bx = a_ql4.box()*a_qr4.box();
  for (int dir = 0; dir < DIM; dir++)
    {
      if (dir!=a_dir) bx = bx
                        *a_ql2.box().grow(dir,-1)*a_qr2.box().grow(dir,-1);
    }
  BoxData<T,C,MEM,D,E> quot(bx);
  for (int ee = 0; ee < E;ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
        {
          auto quotcomp = slice(quot,cc,dd,ee);
          auto qr4comp = slice(a_qr4,cc,dd,ee);
          auto qr2comp = slice(a_qr2,cc,dd,ee);
          Operator::_faceQuotient(quotcomp,a_ql4,qr4comp,a_ql2,qr2comp,a_dir);
        }
  return quot;
}
template<class T, unsigned int CL,unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
BoxData<T,CL,MEM,DR,E>
Operator::_faceMatrixProductAB(
                      const BoxData<T,CL,MEM,DL,E>& a_ql4,
                      const BoxData<T,CR,MEM,DR,E>& a_qr4,
                      const BoxData<T,CL,MEM,DL,E>& a_ql2,
                      const BoxData<T,CR,MEM,DR,E>& a_qr2,
                      int a_dir)
{
  PR_assert(uint(DR) == CL);
  Box bx = a_ql4.box()&a_qr4.box();
  for (int dir = 0; dir < DIM; dir++)
    {
      if (dir!=a_dir) bx = bx
                        *a_ql2.box().grow(dir,-1)*a_qr2.box().grow(dir,-1);
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
template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
BoxData<T,DL,MEM,DR,E>
Operator::_faceMatrixProductATB(
                               const BoxData<T,CL,MEM,DL,E>& a_ql4,
                               const BoxData<T,CR,MEM,DR,E>& a_qr4,
                               const BoxData<T,CL,MEM,DL,E>& a_ql2,
                               const BoxData<T,CR,MEM,DR,E>& a_qr2,
                               int a_dir)
{
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
