
template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
Box BoundaryCondition<T,C,MEM,D,E>::Boundary(Box a_domainBox, Face a_face)
{
    Point dir = Point::Basis(a_face.dir, a_face.sign());
    switch (a_face.side)
    {
        case Side::Lo: return a_domainBox.edge(dir, 1);
        case Side::Hi: return a_domainBox.adjacent(dir, 1);
    }
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
BCStencil<T> BoundaryCondition<T,C,MEM,D,E>::DirichletStencil(T a_dx, Face a_face)
{
    Array<T,5> coefs;
    coefs.fill(1.0/(72.0*a_dx));
    coefs[4] *= -300.0; //multiplies BC
    coefs[0] *= 415.0;
    coefs[1] *= -161.0;
    coefs[2] *= 55.0;
    coefs[3] *= -9.0;

    int sign = a_face.sign();
    BCStencil<T> ret;
    for (int ii = 0; ii < 4; ii++)
    {
        int shift = -sign*(sign == -1 ? ii : ii+1);
        ret.stencil += coefs[ii]*Shift::Basis(a_face.dir, shift);
    }
    ret.coef = coefs[4];
    return ret;
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoundaryCondition<T,C,MEM,D,E>::Dirichlet(
        BoxData<T,C,MEM,D,E>& a_flux,
        const BoxData<T,C,MEM,D,E>& a_state,
        Box a_domainBox,
        T a_dx,
        T a_value,
        Face a_face)
{
    Box bound = Boundary(a_domainBox, a_face);
    bound &= a_flux.box();
    if (bound.empty()) { return; }
    auto S = DirichletStencil(a_dx, a_face);
   
    PROTO_ASSERT(a_state.box().contains(S.stencil.domain(bound)),
            "BoundaryCondition::Dirichlet | Error: Not enough state data to apply BC");

    BoxData<T,C,MEM,D,E> fluxTemp = S.stencil(a_state, bound);
    fluxTemp += (S.coef * a_value);
    fluxTemp *= (-a_face.sign());
    fluxTemp.copyTo(a_flux);
}

template<class T, unsigned int C, MemType MEM, unsigned int D, unsigned int E>
void BoundaryCondition<T,C,MEM,D,E>::Dirichlet(
        BoxData<T,C,MEM,D,E>& a_flux,
        const BoxData<T,C,MEM,D,E>& a_state,
        Box a_domainBox,
        const BoxData<T,C,MEM,D,E>& a_values,
        T a_dx,
        Face a_face)
{
    Box bound = Boundary(a_domainBox, a_face);
    bound &= a_flux.box();
    if (bound.empty()) { return; }
    auto S = DirichletStencil(a_dx, a_face);
   
    PROTO_ASSERT(a_state.box().contains(S.stencil.domain(bound)),
            "BoundaryCondition::Dirichlet | Error: Not enough state data to apply BC");
    PROTO_ASSERT(a_values.box().contains(bound),
            "BoundaryCondition::Dirichlet | Error: Insufficient prescribed BC data");
    BoxData<T,C,MEM,D,E> fluxTemp = S.stencil(a_state, bound);
    Stencil<T> I = 1.0*Shift::Zeros();
    fluxTemp += I(a_values, S.coef);
    fluxTemp *= (-a_face.sign());
    fluxTemp.copyTo(a_flux);
}
