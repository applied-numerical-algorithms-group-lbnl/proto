
template<typename T>
InterpStencil<T>::InterpStencil(Point a_refRatio)
{
    define(a_refRatio);
}

template<typename T>
void InterpStencil<T>::define(Point a_refRatio)
{
    PROTO_ASSERT(a_refRatio.min() > 0, 
        "InterpStencil::define | Error: non-positive refinement ratio");
    m_ratio = a_refRatio;
    for (int dir = 0; dir < DIM; dir++)
    {
        m_interp[dir].define(dir, a_refRatio[dir]); 
    }
}

template<typename T>
Box InterpStencil<T>::domain(Box a_range) const
{
    Point L = Point::Zeros();
    Point H = Point::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    for (int shift = 0; shift < m_interp[dir].ratio(); shift++)
    {
        Box di = m_interp[dir].get(shift).domain(a_range);
        L[dir] = std::min(di.low()[dir], L[dir]);
        H[dir] = std::max(di.high()[dir], H[dir]);
    }
    return Box(L,H);
}

template<typename T>
Box InterpStencil<T>::range(Box a_domain) const
{
    Point L, H;
    for (int dir = 0; dir < DIM; dir++)
    for (int shift = 0; shift < m_interp[dir].ratio(); shift++)
    {
        Box di = m_interp[dir].get(shift).range(a_domain);
        if (shift == 0)
        {
            L[dir] = di.low()[dir];
            H[dir] = di.high()[dir];
        } else {
            L[dir] = std::min(di.low()[dir], L[dir]);
            H[dir] = std::max(di.high()[dir], H[dir]);
        }
    }
    return Box(L,H);
}

template<typename T>
Point InterpStencil<T>::ghost() const
{
    Point span_0 = m_interp[0].ghost();
    for (int jj = 1; jj < DIM; jj++)
    {
        Point span = m_interp[jj].ghost();
        for (int dir = 0; dir < DIM; dir++)
        {
            span_0[dir] = max(span[dir], span_0[dir]);
        }
    }
    return span_0;
}

template<typename T>
Box InterpStencil<T>::span() const
{
    Box b = m_interp[0].span();
    for (int dir = 1; dir < DIM; dir++)
    {
        b = b & m_interp[dir].span().low();
        b = b & m_interp[dir].span().high();
    }
    return b;
}


template<typename T>
template<unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void InterpStencil<T>::apply(
    BoxData<T, C, MEM, D, E>& a_output,
    const BoxData<T, C, MEM, D, E>& a_input,
    Box                       a_range,
    bool                      a_overwrite,
    T                         a_scale) const
{
    PR_TIMERS("InterpStencil::apply");
    Box rangeBox = a_range;
    if (rangeBox.empty())
    {
        rangeBox = range(a_input.box()) & a_output.box();
    }

    Array<BoxData<T, C, MEM>, DIM> stage;
    stage[0] = m_interp[0](a_input);
    for (int dir = 0; dir < DIM-1; dir++)
    {
        stage[dir+1] = m_interp[dir+1](stage[dir]);
    }
    Stencil<T> S = a_scale*Shift::Zeros();
    if (a_overwrite)
    {
        a_output |= S(stage[DIM-1], rangeBox);
    } else {
        a_output += S(stage[DIM-1], rangeBox);
    }
}

template<typename T>
template<unsigned int C, MemType MEM,
    unsigned int D, unsigned int E>
LazyInterpStencil<T, C, MEM, D, E>
InterpStencil<T>::operator()(
        const BoxData<T, C, MEM, D, E>& a_input,
        T a_scale) const
{
    return LazyInterpStencil<T,C,MEM,D,E>(
            *this, a_input, Box(), a_scale);
}

template<typename T>
template<unsigned int C, MemType MEM,
    unsigned int D, unsigned int E>
LazyInterpStencil<T, C, MEM, D, E>
InterpStencil<T>::operator()(
        const BoxData<T, C, MEM, D, E>& a_input,
        Box a_box,
        T a_scale) const
{
#ifndef NDEBUG
    Box b = this->range(a_input.box());
    PROTO_ASSERT(b.containsBox(a_box),
            "Stencil::operator() | Error: User provided Box is larger than the maximum range.");
#endif
    return LazyInterpStencil<T,C,MEM,D,E>(
            *this, a_input, a_box, a_scale);
}

template<typename T, unsigned int C, MemType MEM,
    unsigned int D, unsigned int E>
LazyInterpStencil<T,C,MEM,D,E>::LazyInterpStencil(
        const InterpStencil<T>& a_stencil,
        const BoxData<T, C, MEM, D, E>& a_src,
        Box a_range, T a_scale)
{
    op = &a_stencil;
    range = a_range;
    scale = a_scale;
    srcData = &a_src;
}

template<typename T, unsigned int C, MemType MEM,
    unsigned int D, unsigned int E>
void LazyInterpStencil<T,C,MEM,D,E>::apply(
    BoxData<T,C,MEM,D,E>& a_dst, bool a_overwrite) const
{
    op->apply(a_dst, *srcData, range, a_overwrite, scale);  
}

template<typename T, unsigned int C, MemType MEM,
    unsigned int D, unsigned int E>
Box LazyInterpStencil<T,C,MEM,D,E>::inferredRange() const
{
    return op->range(srcData->box());
}


