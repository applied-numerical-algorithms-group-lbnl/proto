
template<typename T>
InterpStencil<T>::InterpStencil(Point a_refRatio)
{
    define(a_refRatio);
}

template<typename T>
void InterpStencil<T>::define(Point a_refRatio)
{
    for (int dir = 0; dir < DIM; dir++)
    {
        m_interp[dir].define(dir, a_refRatio[dir]); 
    }
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
template<unsigned int C, MemType MEM, unsigned char D, unsigned char E>
void InterpStencil<T>::apply(
    BoxData<T, C, MEM, D, E>& a_output,
    BoxData<T, C, MEM, D, E>& a_input,
    T                   a_scale)
{
  PR_TIMERS("InterpStencil::apply");
    Array<BoxData<T, C, MEM>, DIM+1> stage;
    stage[1] = m_interp[0](a_input);
    for (int dir = 1; dir < DIM; dir++)
    {
        stage[dir+1] = m_interp[dir](stage[dir]);
    }
    Stencil<T> idOp = 1.0*Shift::Zeros();
    a_output |= idOp(stage[DIM]);
}
