template<class T>
void InterpStencil1D<T>::define(unsigned int a_dir, unsigned int a_ratio)
{
    PROTO_ASSERT(a_ratio > 0,
            "InterpStencil1D(uint, uint) | Error: refinement ratio must be positive");
    PROTO_ASSERT(a_dir < DIM,
            "InterpStencil1D(uint, uint) | Error: invalid coordinate %u", m_dir);
    m_r = a_ratio;
    m_dir = a_dir;
    m_s.clear();
    m_s.resize(a_ratio);
}

template<class T>
const Stencil<T>& InterpStencil1D<T>::get(unsigned int a_shift) const
{
    PROTO_ASSERT(a_shift < m_r,
            "InterpStencil1D::get | Error: invalid shift %u", a_shift);
    return m_s[a_shift];
}

template<class T>
void InterpStencil1D<T>::set(Stencil<T> a_S, unsigned int a_shift)
{
    PROTO_ASSERT(a_shift < m_r,
            "InterpStencil1D::get | Error: invalid shift %u", a_shift);
    Point ratio = Point::Ones();
    ratio[m_dir] = m_r;

    a_S.destRatio() = ratio;
    a_S.destShift() = Point::Basis(m_dir, a_shift);
    m_s[a_shift] = a_S;
}

template<class T>
Box InterpStencil1D<T>::span() const
{
    Box span;
    for (auto& S : m_s)
    {
        span = span & S.span().low();
        span = span & S.span().high();
    }
    return span;
}

template<class T>
Point InterpStencil1D<T>::ghost() const
{
    Box spanBox = span();
    int v[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        v[ii] = std::max(std::abs(spanBox.low()[ii]), std::abs(spanBox.high()[ii]));
    }
    return Point(v);
}

template<class T>
void InterpStencil1D<T>::print() const
{
    pr_out() << "InterpStencil1D | r: " << m_r << " | dir: " << m_dir << " | span: " << span() << std::endl;
    for (auto& S : m_s) { m_s.print(); }
}

template<class T>
template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
LazyStencil<T,C,MEMTYPE,D,E> InterpStencil1D<T>::operator()(
        const BoxData<T,C,MEMTYPE,D,E>& a_src,
        Box                             a_box,
        T                               a_scale) const
{
    LazyStencil<T,C,MEMTYPE,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();
    for (int ii = 0; ii < m_s.size(); ii++)
    {
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(&(m_s[ii]));
        ret.m_src.push_back(const_cast<BoxData<T,C,MEMTYPE,D,E>*>(&a_src));
        ret.m_box.push_back(a_box & b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}

template<class T>
template<unsigned int C, MemType MEMTYPE, unsigned char D, unsigned char E>
inline LazyStencil<T,C,MEMTYPE,D,E> InterpStencil1D<T>::operator()(
        const BoxData<T,C,MEMTYPE,D,E>& a_src,
		T                               a_scale) const
{
    LazyStencil<T,C,MEMTYPE,D,E> ret;
    if (m_s.size() <= 0){return ret;}
    Box bsrc = a_src.box();

    for (int ii = 0; ii < m_s.size(); ii++)
    {
        //compute range for building the output BoxData, if necessary
        if (ii == 0)
        {
            ret.m_range = m_s[ii].range(bsrc); //initialize
        } else {
            ret.m_range = ret.m_range & m_s[ii].range(bsrc).low();
            ret.m_range = ret.m_range & m_s[ii].range(bsrc).high();
        }

        //compute the largest possible computation box (different from range) 
        Box b = Box(bsrc.low() - m_s[ii].span().low(), bsrc.high() - m_s[ii].span().high());
        b = b.taperCoarsen(m_s[ii].srcRatio());
        ret.m_stencil.push_back(&(m_s[ii]));
        ret.m_src.push_back(const_cast<BoxData<T,C,MEMTYPE,D,E>*>(&a_src));
        ret.m_box.push_back(b);
        ret.m_scale.push_back(a_scale);
    }
    return ret;
}
