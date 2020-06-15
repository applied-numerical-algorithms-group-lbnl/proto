using namespace Proto;
/*
template<unsigned int C>
void
Operator::faceDivergence(PatchState<C,0>& a_div, const PatchState<0,C>& a_flux, double a_dx)
{
    auto div_i = a_div.cellData();
    for (int ii = 0; ii < DIM; ii++)
    {
        auto flx_i = a_flux.fluxData(ii);
        auto DIV = 1.0*Shift::Basis(ii,1) - 1.0*Shift::Zeros();
        if (ii == 0)
        {
            div_i |= DIV(flx_i, 1.0/a_dx);
        } else {
            div_i += DIV(flx_i, 1.0/a_dx);
        }
    }
}
*/
template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::convolve(BoxData<T,C,D,E>& a_avg,
        const BoxData<T,C,D,E>& a_ctr,
        const BoxData<T,C,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_avg |= Stencil<T>::Laplacian()(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::convolve(const BoxData<T,C,D,E>& a_ctr,
        const BoxData<T,C,D,E>& a_2nd)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,D,E> a_avg(a_ctr.box() & a_2nd.box().grow(-1));
    convolve(a_avg, a_ctr, a_2nd);
    return a_avg;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::deconvolve(BoxData<T,C,D,E>& a_ctr,
        const BoxData<T,C,D,E>& a_avg,
        const BoxData<T,C,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    a_ctr |= Stencil<T>::Laplacian()(a_2nd, -1.0/24);
    a_ctr += a_avg;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::deconvolve(const BoxData<T,C,D,E>& a_avg,
        const BoxData<T,C,D,E>& a_2nd)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1)),
            "error Operator::convolve | centered data defined on too small a box.");

    BoxData<T,C,D,E> a_ctr(a_avg.box() & a_2nd.box().grow(-1));
    deconvolve(a_ctr, a_avg, a_2nd);
    return a_ctr;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::convolveFace(BoxData<T,C,D,E>& a_avg,
        const BoxData<T,C,D,E>& a_ctr,
        const BoxData<T,C,D,E>& a_2nd,
        int a_dir)
{
    PROTO_ASSERT(a_ctr.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_avg |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, 1.0/24.0);
    a_avg += a_ctr;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::convolveFace(const BoxData<T,C,D,E>& a_ctr,
        const BoxData<T,C,D,E>& a_2nd,
        int a_dir)
{
    auto b = a_2nd.box().grow(-1).grow(a_dir, 1);
    b &= a_ctr.box();
    BoxData<T,C,D,E> avg(b);
    convolveFace(avg, a_ctr, a_2nd, a_dir);
    return avg;
}


template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::deconvolveFace(BoxData<T,C,D,E>& a_ctr,
        const BoxData<T,C,D,E>& a_avg,
        const BoxData<T,C,D,E>& a_2nd,
        int a_dir)
{
    PROTO_ASSERT(a_avg.box().contains(a_2nd.box().grow(-1).grow(a_dir,1)),
    "Error in Operator::convolveFace | Insufficient source data.");
    a_ctr |= Stencil<T>::LaplacianFace(a_dir)(a_2nd, -1.0/24.0);
    a_ctr += a_avg;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::deconvolveFace(const BoxData<T,C,D,E>& a_avg,
        const BoxData<T,C,D,E>& a_2nd,
        int a_dir)
{
    auto b = a_2nd.box().grow(-1).grow(a_dir, 1);
    b &= a_avg.box();
    BoxData<T,C,D,E> ctr(b);
    deconvolveFace(ctr, a_avg, a_2nd, a_dir);
    return ctr;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::cellProduct(
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd)
{
    auto b = (a_p2nd.box().grow(-1) & a_q2nd.box().grow(-1));
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,D,E> pq(b);
    cellProduct(pq, a_pavg, a_qavg, a_p2nd, a_q2nd);
    return pq;
}
template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::cellProduct(BoxData<T,C,D,E>& a_product,
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellProduct | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellProduct | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellProduct | incompatible patch sizes.");
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,D,E> Dq = Diff(a_q2nd);

        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_dp, Var<T,C,D,E>& v_dq, int v_dir)
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
        [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_p, Var<T,C,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += v_p(cc,dd,ee)*v_q(cc,dd,ee);
            }
        }, a_product, a_pavg, a_qavg);
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::faceProduct(
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd,
        int a_dir)
{
    Proto::Box b = (a_p2nd.box() & a_q2nd.box()).grow(-1).grow(a_dir, 1);
    
    BoxData<T,C,D,E> out(b);
    faceProduct(out, a_pavg, a_qavg, a_p2nd, a_q2nd, a_dir);
    return out;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::faceProduct(BoxData<T,C,D,E>& a_dest,
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd,
        int a_dir)
{
    bool do_replace = true;
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){continue;}
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_dp, Var<T,C,D,E>& v_dq, bool v_replace)
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
        [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_p, Var<T,C,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += v_p(cc,dd,ee)*v_q(cc,dd,ee);
            }
        }, a_dest, a_pavg, a_qavg);
}


template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::cellQuotient(
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd)
{
    auto b = (a_p2nd.box().grow(-1) & a_q2nd.box().grow(-1));
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,D,E> quotient(b);
    cellQuotient(quotient, a_pavg, a_qavg, a_p2nd, a_q2nd);
    return quotient;

}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::cellQuotient(BoxData<T,C,D,E>& a_quotient,
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    //BoxData<T,C,D,E> out(a_p2nd.box().grow(-1));
    //out.setVal(0);
    for (int ii = 0; ii < DIM; ii++)
    {
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_dp, Var<T,C,D,E>& v_dq, Var<T,C,D,E>& v_p, Var<T,C,D,E>& v_q, int v_dir)
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
        [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_p, Var<T,C,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::faceQuotient(
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd,
        int a_dir)
{
    Proto::Box b = (a_p2nd.box() & a_q2nd.box()).grow(-1).grow(a_dir, 1);
    b &= (a_pavg.box() & a_qavg.box());
    BoxData<T,C,D,E> out(b);
    faceQuotient(out, a_pavg, a_qavg, a_p2nd, a_q2nd, a_dir);
    return out;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::faceQuotient(BoxData<T,C,D,E>& a_quotient,
        const BoxData<T,C,D,E>& a_pavg,
        const BoxData<T,C,D,E>& a_qavg,
        const BoxData<T,C,D,E>& a_p2nd,
        const BoxData<T,C,D,E>& a_q2nd,
        int a_dir)
{
    PROTO_ASSERT(a_p2nd.box() == a_q2nd.box(), "Operator::cellQuotient | p2nd and q2nd should have the same box.");
    PROTO_ASSERT(a_pavg.box().contains(a_p2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    PROTO_ASSERT(a_qavg.box().contains(a_q2nd.box().grow(-1)), "Operator::cellQuotient | incompatible patch sizes.");
    
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){continue;}
        auto Diff = 0.5*Shift::Basis(ii,1) - 0.5*Shift::Basis(ii,-1);
        BoxData<T,C,D,E> Dp = Diff(a_p2nd);
        BoxData<T,C,D,E> Dq = Diff(a_q2nd);
        forallInPlace(
            [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_dp, Var<T,C,D,E>& v_dq, Var<T,C,D,E>& v_p, Var<T,C,D,E>& v_q, int v_dir)
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
        [=] PROTO_LAMBDA (Var<T,C,D,E>& v_out, Var<T,C,D,E>& v_p, Var<T,C,D,E>& v_q)
        {
            for (int cc = 0; cc < C; cc++)
            for (int dd = 0; dd < D; dd++)
            for (int ee = 0; ee < E; ee++)
            {
                v_out(cc,dd,ee) += (v_p(cc,dd,ee) / v_q(cc,dd,ee));
            }
        }, a_quotient, a_pavg, a_qavg);
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::cellToFace(const BoxData<T,C,D,E>& a_cell, int a_dir)
{
    auto S = Proto::Stencil<double>::CellToEdge(a_dir);
    BoxData<T,C,D,E> face = S(a_cell);
    return face;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::cellToFace(BoxData<T,C,D,E>& a_face,
    const BoxData<T,C,D,E>& a_cell,
    int a_dir)
{
    auto S = Proto::Stencil<double>::CellToEdge(a_dir);
    a_face |= S(a_cell);
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
BoxData<T,C,D,E>
Operator::diffCellToFace(const BoxData<T,C,D,E>& a_data, int a_dir)
{
    auto Diff = Proto::Stencil<double>::DiffCellToEdge(a_dir);
    BoxData<T,C,D,E> deriv = Diff(a_data);
    return deriv;
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
void
Operator::diffCellToFace(BoxData<T,C,D,E>& a_deriv,
    const BoxData<T,C,D,E>& a_data,
    int a_dir)
{
    auto Diff = Proto::Stencil<double>::DiffCellToEdge(a_dir);
    a_deriv |= Diff(a_data);
}



