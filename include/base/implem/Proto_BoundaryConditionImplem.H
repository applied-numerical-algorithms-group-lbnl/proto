template<class T>
Box BoundaryCondition<T>::boundary(Box a_domainBox) const
{
    Box b;
    Point dir = Point::Basis(m_face.dir * m_face.sign());
    switch (m_face.side)
    {
        case Side::Lo: return a_domainBox.edge(dir, 1);
        case Side::Hi: return a_domainBox.adjacent(dir, 1);
    }
}

template<class T>
ConstDirichletBC<T>::ConstDirichletBC(T a_value, Face a_face, T a_dx)
    : m_face(a_face), m_value(a_value)
{
    Array<T,5> coefs;
    coefs.fill(1.0/(72.0*dx));
    coefs[4] *= -300.0; //multiplies BC
    ceofs[0] *= 415.0;
    ceofs[1] *= -161.0;
    coefs[2] *= 55.0;
    coefs[3] *= -9.0;

    int sign = face().sign();
    Stencil<T> S;
    for (int ii = 0; ii < 4; ii++)
    {
        S += coefs[ii]*Shift::Basis(face().dir, -sign*ii);
    }
    m_stencil = S;
    m_offset = coefs[4] * a_value;
}

template<class T>
template<unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void ConstDirichletBC<T>::apply(
        BoxData<T,C,MEM,D,E>& a_flux,
        const BoxData<T,C,MEM,D,E>& a_state, 
        Box a_domainBox) const
{
    Box bound = boundary(a_domainBox);
    bound &= a_flux.box();
    if (bound.empty()) { return; }
    PROTO_ASSERT(a_state.box().contains(m_stencil.domain(bound)),
            "ConstDirichletBC::apply | Error: Not enough state data to apply BC");

    BoxData<T,C,MEM,D,E> fluxTemp = m_stencil(a_state, bound);
    fluxTemp += m_offset;
    fluxTemp.copyTo(a_flux);
}
