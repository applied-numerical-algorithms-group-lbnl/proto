template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PR_TIME("Operator::convolve");
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_avg |= Stencil<T>::Laplacian()(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
convolve(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr)
{
    _convolve(a_avg, a_ctr, a_ctr);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
_convolve(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd)
{
    PR_TIME("Operator::deconvolve");
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_ctr |= Stencil<T>::Laplacian()(a_2nd, -1.0/24);
    a_ctr += a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
deconvolve(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg)
{
    _deconvolve(a_ctr, a_avg, a_avg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
_deconvolve(
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    PR_TIME("Operator::convolveFace");
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_avg |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
convolveFace(
        BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    _convolveFace(a_avg, a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
convolveFace(
        const BoxData<T,C,MEM,D,E>& a_ctr,
        int a_dir)
{
    return _convolveFace(a_ctr, a_ctr, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        const BoxData<T,C,MEM,D,E>& a_2nd,
        int a_dir)
{
    PR_TIME("Operator::deconvolveFace");
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_ctr |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, -1.0/24.0);
    a_ctr += a_avg;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
deconvolveFace(
        BoxData<T,C,MEM,D,E>& a_ctr,
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    _deconvolveFace(a_ctr, a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
deconvolveFace(
        const BoxData<T,C,MEM,D,E>& a_avg,
        int a_dir)
{
    return _deconvolveFace(a_avg, a_avg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
cellProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellProduct(a_pavg, a_qavg, a_pavg, a_qavg);
}
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
    PR_TIME("Operator::cellProduct");
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);

        forallInPlace(
            [] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_dp, Var<T,C,MEM,D,E>& v_dq, int v_dir)
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
        [] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
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
cellProduct(BoxData<T,C,MEM,D,E>& a_product,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellProduct(a_product, a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
faceProduct(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceProduct(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
    PR_TIME("Operator::faceProduct");
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
faceProduct(BoxData<T,C,MEM,D,E>& a_dest,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceProduct(a_dest, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}


template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
cellQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    return _cellQuotient(a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd)
{
  PR_TIME("Operator::cellQuotient");
  //PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    //BoxData<T,C,MEM,D,E> out(a_p2nd.box().grow(-1));
    //out.setVal(0);
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) = (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);

    BoxData<T,C,MEM,D,E> dqsqFactor(a_quotient.box());
    BoxData<T,C,MEM,D,E> dpdqTerm(a_quotient.box());
    
    dqsqFactor.setVal(1.0);
    dpdqTerm.setToZero();
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
        
        forallInPlace
          ([=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out1,
                             Var<T,C,MEM,D,E>& v_out2,
                             Var<T,C,MEM,D,E>& v_dp,
                             Var<T,C,MEM,D,E>& v_dq,
                             Var<T,C,MEM,D,E>& v_q)
            {
                for (int cc = 0; cc < C; cc++)
                for (int dd = 0; dd < D; dd++)
                for (int ee = 0; ee < E; ee++)
                {                  
                  v_out1(cc,dd,ee) += v_dq(cc,dd,ee)*v_dq(cc,dd,ee)/
                    (v_q(cc,dd,ee)*v_q(cc,dd,ee)*12.0);
                  v_out2(cc,dd,ee) -= v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/
                     (v_q(cc,dd,ee)*v_q(cc,dd,ee)*12.0);
                }
            }, dqsqFactor,dpdqTerm, Dp, Dq, a_qavg);
    }
    a_quotient*=dqsqFactor;
    a_quotient += dpdqTerm;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
cellQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg)
{
    _cellQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
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

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
faceQuotient(
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    return _faceQuotient(a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
_faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        const BoxData<T,C,MEM,D,E>& a_p2nd,
        const BoxData<T,C,MEM,D,E>& a_q2nd,
        int a_dir)
{
  PR_TIME("Operator::faceQuotient");
  //PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    //BoxData<T,C,MEM,D,E> out(a_p2nd.box().grow(-1));
    //out.setVal(0);
    forallInPlace(
        [=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out, Var<T,C,MEM,D,E>& v_p, Var<T,C,MEM,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) = (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);

    BoxData<T,C,MEM,D,E> dqsqFactor(a_quotient.box());
    BoxData<T,C,MEM,D,E> dpdqTerm(a_quotient.box());
    
    dqsqFactor.setVal(1.0);
    dpdqTerm.setToZero();
    for (int ii = 0; ii < DIM; ii++)
    {
      if (ii != a_dir)
        {
          auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
          BoxData<T,C,MEM,D,E> Dp = Diff(a_p2nd);
          BoxData<T,C,MEM,D,E> Dq = Diff(a_q2nd);
          
          forallInPlace
            ([=] PROTO_LAMBDA (Var<T,C,MEM,D,E>& v_out1,
                               Var<T,C,MEM,D,E>& v_out2,
                               Var<T,C,MEM,D,E>& v_dp,
                               Var<T,C,MEM,D,E>& v_dq,
                               Var<T,C,MEM,D,E>& v_q)
             {
               for (int cc = 0; cc < C; cc++)
                 for (int dd = 0; dd < D; dd++)
                   for (int ee = 0; ee < E; ee++)
                     {                  
                       v_out1(cc,dd,ee) += v_dq(cc,dd,ee)*v_dq(cc,dd,ee)/
                         (v_q(cc,dd,ee)*v_q(cc,dd,ee)*12.0);
                       v_out2(cc,dd,ee) -= v_dp(cc,dd,ee)*v_dq(cc,dd,ee)/
                         (v_q(cc,dd,ee)*v_q(cc,dd,ee)*12.0);
                     }
             }, dqsqFactor,dpdqTerm, Dp, Dq, a_qavg);
        }
    }
    a_quotient*=dqsqFactor;
    a_quotient += dpdqTerm;
}

template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void
faceQuotient(BoxData<T,C,MEM,D,E>& a_quotient,
        const BoxData<T,C,MEM,D,E>& a_pavg,
        const BoxData<T,C,MEM,D,E>& a_qavg,
        int a_dir)
{
    _faceQuotient(a_quotient, a_pavg, a_qavg, a_pavg, a_qavg, a_dir);
}

template <class T, MemType MEM>
BoxData<T,DIM,MEM>
cofactor(
        const BoxData<T,DIM,MEM>& a_X,
        int a_dir)
{
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
    int dirperp1 = (a_dir+1)%DIM;
    int dirperp2 = (a_dir+2)%DIM;

    // average X in direction 1
    BoxData<T,DIM,MEM> XAv1 = Stencil<T>::faceToCell(dirperp1,4)(a_X);
    // dX/dx1
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
    auto ret = forall<T,DIM,MEM>([]PROTO_LAMBDA(Var<T,DIM,MEM>& a_ret, Var<T,DIM,MEM>& a_d21,Var<T,DIM,MEM>& a_d12)
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

template <class T, MemType MEM>
BoxData<T,1,MEM>
jacobian(
        const BoxData<T,DIM,MEM>& a_X,
        const std::array<BoxData<T,DIM,MEM>, DIM>& a_NT)
{
    PR_TIME("Operator::jacobian");
    std::array<BoxData<T,1,MEM>,DIM> dfdx;
    for (int dir = 0; dir < DIM; dir++)
    {
        BoxData<T,DIM,MEM> xFace4 = Stencil<T>::cornersToFaces(dir,4)(a_X);
        BoxData<T,DIM,MEM> xFace2 = Stencil<T>::cornersToFaces(dir,2)(a_X);
        auto flux = Operator::_faceMatrixProductATB(a_NT[dir],xFace4,a_NT[dir],xFace2,dir);
        dfdx[dir] = Stencil<T>::FluxDivergence(dir)(flux);
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
  std::array<BoxData<T,1,MEM>,DIM> qlcomp2,qrcomp2,qlcomp4,qrcomp4;  
  for (int comp = 0; comp < DIM; comp++)
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
    forall<T,DIM,MEM>([]PROTO_LAMBDA(Var<T,DIM,MEM>& a_retval,
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
template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
BoxData<T,C,MEM,D,E>
_product1D(
                    const BoxData<T,C,MEM,D,E>& a_ql4,
                    const BoxData<T,C,MEM,D,E>& a_qr4,
                    const BoxData<T,C,MEM,D,E>& a_ql2,
                    const BoxData<T,C,MEM,D,E>& a_qr2,
                    int a_dir)
{
  PR_TIME("Operator::product1D");
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
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E>
_faceTensorQuotient(
                              const BoxData<T,C,MEM,D,E>& a_ql4,
                              const BoxData<T,1,MEM,1,1>& a_qr4,
                              const BoxData<T,C,MEM,D,E>& a_ql2,
                              const BoxData<T,1,MEM,1,1>& a_qr2,        
                              int a_dir)
{
  PR_TIME("Operator::faceTensorQuotient");
  Box bx = a_ql4.box()&a_qr4.box();
  for (int dir = 0; dir < DIM; dir++)
    {
      if (dir!=a_dir) bx = bx
                        &a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
    }
  BoxData<T,C,MEM,D,E> quot(bx);
  for (int ee = 0; ee < E;ee++)
    for (int dd = 0; dd < D; dd++)
      for (int cc = 0; cc < C; cc++)
        {
          auto quotcomp = slice(quot,cc,dd,ee);
          auto ql4comp = slice(a_ql4,cc,dd,ee);
          auto ql2comp = slice(a_ql2,cc,dd,ee);
          Operator::_faceQuotient(quotcomp,ql4comp,a_qr4,ql2comp,a_qr2,a_dir);
        }
  return quot;
}
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
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
      if (dir!=a_dir) bx = bx
                        &a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
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
template<class T, unsigned int CL,unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
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
  Box bx = a_ql4.box()&a_qr4.box();
  for (int dir = 0; dir < DIM; dir++)
    {
      if (dir!=a_dir) bx = bx
                        &a_ql2.box().grow(dir,-1)&a_qr2.box().grow(dir,-1);
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
template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
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
template<class T, unsigned int CL,unsigned char DL,
           unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
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
template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
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
template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
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
template<class T, MemType MEM>
BoxData<T,DIM,MEM,DIM>
_cofactorMatrix(
        const array<BoxData<T,DIM,MEM>, DIM>& a_NTFaces,
        int a_dirNorm)                         
{
  PR_TIME("Operator::cofactorMatrix");
  array<BoxData<T,DIM,MEM>, DIM> NTAveToFace;
  for (int dir = 0; dir < DIM; dir++)
    {
      if (dir != a_dirNorm)
        {
          NTAveToFace[dir] =
            (Stencil<T>::faceToCell(dir,4)*
             Stencil<T>::CellToFace(a_dirNorm))(a_NTFaces[dir]);
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
template<class T,unsigned int C, MemType MEM>
BoxData<T,C,MEM,DIM>
_faceGradient(
              const BoxData<T,C,MEM>& a_phi,
              int a_dirNorm,
              int a_order)
{
    PR_TIME("Operator::faceGradient");
  // This is an operator, not a stencil, since there are different stencils used to
  // Compute different components of the gradient.
  PROTO_ASSERT((a_order == 2) || (a_order == 4),
               "_faceGradient defined only for orders = 2 and 4");
  array<Stencil<T>,DIM> faceDerivs; 
  for (int dir = 0;dir < DIM;dir++)
    {
      if (dir == a_dirNorm)
        {
          faceDerivs[dir] = Stencil<T>::DiffCellToFace(a_dirNorm,Side::Lo,a_order);
        }
      else
        {
          faceDerivs[dir] =           
            Stencil<T>::CellToFace(dir,Side::Lo,a_order)*
            Stencil<T>::FluxDivergence(dir)*
            Stencil<T>::CellToFace(a_dirNorm,Side::Lo,a_order);
        }
    }
  array<BoxData<T,C,MEM>,DIM> gradComps;
  for (int dir = 0; dir < DIM; dir++)
    {
          gradComps[dir] = faceDerivs[dir](a_phi);
    }
#if DIM==2
  auto retval =
    forall<T,C,MEM,DIM,1>
    ([ ] PROTO_LAMBDA(
                      Var<T,C,MEM,DIM>& a_retval,
                      const Var<T,C,MEM>& a_dphidx0,
                      const Var<T,C,MEM>& a_dphidx1)
     {
       for (int cc = 0;cc < C;cc++)
         {
           a_retval(cc,0) = a_dphidx0(cc);
           a_retval(cc,1) = a_dphidx1(cc);
         }
     },
     gradComps[0],gradComps[1]); 
#elif DIM==3
  auto retval =
    forall<T,C,MEM,DIM,1>
    ([ ] PROTO_LAMBDA(
                      Var<T,C,MEM,DIM>& a_retval,
                      const Var<T,C,MEM>& a_dphidx0,
                      const Var<T,C,MEM>& a_dphidx1,
                      const Var<T,C,MEM>& a_dphidx2)
     {
       for (int cc = 0; cc < C;cc++)
         {
           a_retval(cc,0) = a_dphidx0(cc);
           a_retval(cc,1) = a_dphidx1(cc);
           a_retval(cc,2) = a_dphidx2(cc);
         }
     },
     gradComps[0],gradComps[1],gradComps[2]);
#endif
  return retval;
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
