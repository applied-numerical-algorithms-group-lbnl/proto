#ifndef PRroto_stencildefs
#define PRroto_stencildefs
#include "Proto_Stencil.H"

//========================================================================
// Stencil Definitions
//======================
template<typename T>
Stencil<T> Stencil<T>::Laplacian()
{
    Stencil<T> ret = ((T)(-2*DIM))*Shift::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
      ret += ((T)(1))*Shift::Basis(dir, 1);
      ret += ((T)(1))*Shift::Basis(dir,-1);
    }
    return ret;
}
#if DIM == 2
template<typename T>
Stencil<T> Stencil<T>::Laplacian_9()
{
  T factor = -20./6.;
    Stencil<T> ret = factor*Shift::Zeros();
    Box K = Box::Cube(3).shift(Point::Ones(-1));
    T a = 1.0/6.0;
    T b = 2.0/3.0;
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Point p = *iter;
        int pointSum = abs(p[0]) + abs(p[1]);
        if (pointSum == 1)
        {
            ret += b*Shift(p); 
        } else if (pointSum == 2)
        {
            ret += a*Shift(p); 
        }
    }
    return ret;
}
#elif DIM == 3
template<typename T>
Stencil<T> Stencil<T>::Laplacian_19()
{
    Stencil<T> ret = ((T)(-24.0))*Shift::Zeros();
    Box K = Box::Cube(3).shift(Point::Ones(-1));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        int sum = (*iter).absSum();
        if (sum == 1) {
            ret += ((T)(2.0))*Shift(*iter);
        } else if (sum == 2) {
            ret += ((T)(1.0))*Shift(*iter);
        }
    }
    ret *= (1.0/6.0);
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::Laplacian_27()
{
    Stencil<T> ret = ((T)(-64.0/15.0))*Shift::Zeros();
    Box K = Box::Cube(3).shift(Point::Ones(-1));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Point p = *iter;
        int pointSum = p.absSum();
        if (pointSum == 1)
        {
            ret += ((T)(7.0/15.0))*Shift(p);
        } else if (pointSum == 2) {
            ret += ((T)(1.0/10.0))*Shift(p);
        } else if (pointSum == 3) {
            ret += ((T)(1.0/30.0))*Shift(p);
        }
    }
    return ret;
}
#endif
template<typename T>
Stencil<T> Stencil<T>::LaplacianFace(int a_dir, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::Laplacian(int a_order, int a_dir) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
    switch (a_order)
    {
        case 2:
        {
            ret = ((T)(-2*(DIM-1)))*Shift(Point::Zeros());
            for (int dir = 0; dir < DIM; dir++)
            {
                if (a_dir != dir)
                {
                    ret += ((T)(1))*Shift(Point::Basis(dir));
                    ret += ((T)(1))*Shift(Point::Basis(dir)*(-1));
                }
            }
            break;
        }
        default:
        {
            std::cout << "Stencil::Laplacian(int a_order, int a_dir) invalid for a_order = ";
            std::cout << a_order << ". This configuration is not yet implemented as a built-in Stencil." << std::endl;
            std::abort();
            break;
        }
    }
    return ret;
}

template<typename T>
Stencil<T>
Stencil<T>::CellToEdge(int a_dir, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::CellToEdgeL(int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
//    Point edir = Point::Basis(a_dir);
    switch (a_order)
    {
        case 4:
        {
            T coefs[4] = {-1,7,7,-1};
            ret = 
                coefs[0]*Shift::Basis(a_dir, -2) + coefs[3]*Shift::Basis(a_dir, 1) +
                coefs[1]*Shift::Basis(a_dir, -1)  + coefs[2]*Shift::Zeros();
            ret *= (T)(1.0/12.0);
            break;
        }
        case 5:
        {
            T coefs[5] = {-3,27,47,-13,2};
            ret =  
                coefs[0]*Shift::Basis(a_dir, -2) + coefs[4]*Shift::Basis(a_dir, 2) +
                coefs[1]*Shift::Basis(a_dir, -1) + coefs[3]*Shift::Basis(a_dir, 1) +
                coefs[2]*Shift::Zeros();
            ret *= (T)(1/60.0);
            break;
        }
        default:
        {
            std::abort();
            break;
        }
    }
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::DiffCellToEdge(int a_dir, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::DiffCellToEdge(int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
    switch (a_order)
    {
        case 4:
        {
            T coefs[4] = {1, -15, 15, -1};
            ret = coefs[0]*Shift::Basis(a_dir, -2) +
                  coefs[1]*Shift::Basis(a_dir, -1) +
                  coefs[2]*Shift::Zeros() +
                  coefs[3]*Shift::Basis(a_dir, 1);
            ret *= (T)(1.0/12.0);
            break;
        }
        default:
        {
            std::cout << "Stencil::DiffCellToEdge(int a_dir, int a_order) invalid for a_order = ";
            std::cout << a_order << ". This configuration is not yet implemented as a built-in Stencil." << std::endl;
            std::abort();
            break;
        }
    }
}

template<typename T>
Stencil<T> Stencil<T>::CellToEdgeL(int a_dir, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::CellToEdgeL(int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    Stencil<T> ret;
    //PROTO_ASSERT(false, "Call to Stencil::CellToEdgeL(...) | This method is deprecated. Please use Stencil::CellToEdge");
    
    switch (a_order)
    {
        case 5:
        {
            Point edir = Point::Basis(a_dir);
            T coefs[5] = {2,-13,47,27,-3};
            ret =  
                coefs[0]*Shift(edir*(-3)) + coefs[4]*Shift(edir*(1)) +
                coefs[1]*Shift(edir*(-2)) + coefs[3]*Shift(Point::Zeros()) +
                coefs[2]*Shift(edir*(-1));
            ret *= (T)(1/60.0);
            break;
        }
        default:
        {
            std::cout << "Stencil::CellToEdgeL(int a_dir, int a_order) invalid for a_order = ";
            std::cout << a_order << ". This configuration is not yet implemented as a built-in Stencil." << std::endl;
            std::abort();
            break;
        }
    }
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::CellToEdgeH(int a_dir, int a_order)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::CellToEdgeH(int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    //PROTO_ASSERT(false, "Call to Stencil::CellToEdgeH(...) | This method is deprecated. Please use Stencil::CellToEdge");
    Stencil<T> ret;
    switch (a_order)
    {
        case 5:
        {
            Point edir = Point::Basis(a_dir);
            T coefs[5] = {-3,27,47,-13,2};
            ret =  
                coefs[0]*Shift(edir*(-2)) + coefs[4]*Shift(edir*(2)) +
                coefs[1]*Shift(edir*(-1)) + coefs[3]*Shift(edir*(1)) +
                coefs[2]*Shift(Point::Zeros());
            ret *= (T)(1/60.0);
            break;
        }
        default:
        {
            std::cout << "Stencil::CellToEdgeL(int a_dir, int a_order) invalid for a_order = ";
            std::cout << a_order << ". This configuration is not yet implemented as a built-in Stencil." << std::endl;
            std::abort();
            break;
        }
    }
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::EdgeToCell(int a_dir, int a_order)
{    
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::EdgeToCell(int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    PROTO_ASSERT(false, "Call to Stencil::EdgeToCell(...) | This method is deprecated (and was never quite right). Please use Stencil::CellToEdge");
    Stencil<T> ret;
    switch (a_order)
    {
        case 4:
        {
            Point edir = Point::Basis(a_dir);
            T coefs[4] = {-1,7,7,-1};
            ret = 
                coefs[0]*Shift(edir*(-1)) + coefs[3]*Shift(edir*2) +
                coefs[1]*Shift(Point::Zeros()) + coefs[2]*Shift(edir);
            ret *= (T)(1.0/12.0);
            break;
        }
        default:
        {
            std::cout << "No formula found for EdgeToCell of order " << a_order << std::endl;
            break;
        }
    }
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::AvgDown(int a_refRatio)
{
    PROTO_ASSERT(a_refRatio > 0,
    "Stencil::AvgDown(int a_refRatio) invalid for a_refRatio = %i. a_refRatio should be positive.",
    a_refRatio);
    Stencil<T> ret;

    T coef = 1.0/ipow<DIM>(a_refRatio);
    Box kernel = Box::Cube(a_refRatio);
    for (auto iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
        ret += coef*Shift(*iter);
    }
    ret.srcRatio() = Point::Ones(a_refRatio);
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::AvgDownEdge(int a_dir, int a_refRatio)
{
    PROTO_ASSERT(a_refRatio > 0,
    "Stencil::AvgDownEdge(int a_refRatio) invalid for a_refRatio = %i. a_refRatio should be positive.",
    a_refRatio);
    PROTO_ASSERT(DIM > 1, "Error in Stencil::AvgDownEdge | Stencil is invalid for DIM == 1.");
    Stencil<T> ret;

    T coef = 1.0/ipow<DIM-1>(a_refRatio);
    int tuple[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == a_dir){tuple[ii] = 0;}
        else {tuple[ii] = a_refRatio - 1;}
    }
    Point p(tuple); 
    Box kernel(p);
    for (auto iter = kernel.begin(); iter != kernel.end(); ++iter)
    {
        ret += coef*Shift(*iter);
    }
    ret.srcRatio() = Point::Ones(a_refRatio);
    return ret;
}

template<typename T>
Stencil<T> Stencil<T>::FluxDivergence(int a_dir)
{
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::FluxDivergence(int a_dir) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    return ((T)(-1))*Shift(Point::Zeros()) + ((T)(1.0))*Shift(Point::Basis(a_dir));
}

template<typename T>
Stencil<T> Stencil<T>::Derivative(int a_n, int a_dir, int a_order)
{
  static const int MAXORDH = 6;
    PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_dir = %i.\
    a_dir should be in [0,DIM=%i).",a_dir, DIM);
    PROTO_ASSERT((a_n > 0) && (a_n <= 2*MAXORDH),
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_n = %i.\
    a_n must be in [1,%i].",a_n,2*MAXORDH);
    PROTO_ASSERT((a_order > 0) && (a_order <= 2*MAXORDH),
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_order = %i.\
    a_order must be in [1,%i].",a_order,2*MAXORDH);
    PROTO_ASSERT(2*MAXORDH +2 - a_n - a_order >= 0,
    "Stencil::Derivative(int a_n, int a_dir, int a_order) invalid for a_n = %i and a_order = %i.\
    Implementation of requested accuracy not available.\n The maximum accuracy available for a_n = %i is\
    %i.",a_n, a_order, a_n, 2*MAXORDH +2 - (a_n/2)*2);

    if (a_order % 2 != 0){a_order ++;}

    Point edir = Point::Basis(a_dir);
    Stencil<T> factor = ((T)(1.0))*Shift(edir) + ((T)(-2.0))*Shift(Point::Zeros()) + ((T)(1.0))*Shift(edir*(-1.0));
    Stencil<T> deriv;
    Stencil<T> temps[MAXORDH];
    if (a_n % 2 == 0) //even derivative
    {
        T evenMatrixCoefs[MAXORDH][MAXORDH] =
          {{1.,-1./12,1./90,-1./560,1./3150,-1./16632},
           {0,1.,-1./6,7./240,-41./7560,479./453600},
           {0,0,1.,-1./4,13./240,-139./12096},
           {0,0,0,1.,-1./3,31./360},{0,0,0,0,1.,-5./12},
           {0,0,0,0,0,1.}};
        
        int p = (a_n - 2)/2;
        temps[0] = factor;
        for (int ord = 1; ord < MAXORDH; ord++)
        {
            temps[ord] = temps[ord-1]*factor;
        }
        for (int ii = 0; ii < (a_n + a_order - 2)/2; ii++)
        {
            deriv += evenMatrixCoefs[p][ii]*temps[ii];
        }
    } else {//odd derivative
        T oddMatrixCoefs[MAXORDH][MAXORDH] =
          {{1.,-1./6,1./30,-1./140,1./630,-1./2772},
           {0,1.,-1./4,7./120,-41./3024,479./151200},
           {0,0,1.,-1./3,13./144,-139./6048},
           {0,0,0,1.,-5./12,31./240},
           {0,0,0,0,1.,-1./2},
           {0,0,0,0,0,1.}};
        
        int p = (a_n - 1)/2;
        temps[0] = ((T)(0.5))*Shift(edir) + ((T)(-0.5))*Shift(edir*(-1));
        for (int ord = 1; ord < MAXORDH; ord++)
        {
            temps[ord] = temps[ord-1]*factor;
        }
        for (int ii = 0; ii < (a_n + 1 + a_order - 2)/2; ii++)
        {
            deriv += oddMatrixCoefs[p][ii]*temps[ii];
        }
    }
    return deriv;
    
}

//========================================================================
// InterpStencil Definitions
//===========================

template<class T>
InterpStencil<T> InterpStencil<T>::PiecewiseConstant(Point a_ratio)
{
    InterpStencil<T> ret(a_ratio);
    Box K = Box(a_ratio - Point::Ones());
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        Stencil<T>& S = ret(*iter);
        S = 1.0*Shift::Zeros();
        S.destRatio() = a_ratio;
        S.destShift() = (*iter);
    }
    return ret;
}

template<class T>
InterpStencil<T> InterpStencil<T>::PiecewiseLinear(Point a_ratio)
{
    InterpStencil<T> ret(a_ratio);
    Box K = Box(a_ratio - Point::Ones());
    Box U = Box::Cube(2);
    for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
    {
        Point p = *kiter;
        Point q = a_ratio - p;
        Stencil<double>& S = ret(*kiter);
        for (auto uiter = U.begin(); uiter != U.end(); ++ uiter)
        {
            Point e = *uiter;
            double coef = 1.0;
            for (int ii = 0; ii < DIM; ii++)
            {
                if (e[ii] == 0)
                {
                    coef *= q[ii]/((double)a_ratio[ii]);
                } else {
                    coef *= p[ii]/((double)a_ratio[ii]);
                }
            }
            if (coef != 0)
            {
                S += coef*Shift(e);
            }
        }
        S.destRatio() = a_ratio;
        S.destShift() = p;
    }
    return ret;    
}

// helper function for Quadratic InterpStencil
inline double qcoef(int a_ni, int a_r, int a_side)
{
    PROTO_ASSERT(std::abs(a_side) <= 1,
    "qcoef(ni, r, side) invalide for side = %i. Value must be in [-1,1]", a_side);
    double n = a_ni; double r = a_r;
    double coef = 0;
    if (a_side == -1)
    {
        coef = 1.0 + 3.0*n + 3.0*n*n - 3.0*r - 6.0*n*r + 2.0*r*r;
    } 
    else if (a_side == 0)
    {
        coef = -2.0 - 6.0*n - 6.0*n*n + 3.0*r + 6.0*n*r + 5.0*r*r;
    }
    else if (a_side == 1)
    {
        coef = -(-1.0 - 3.0*n - 3.0*n*n + r*r);
    }
    coef /= (6.0*r*r);
    return coef;
}

template<class T>
InterpStencil<T> InterpStencil<T>::Quadratic(int a_ratio)
{
    //can make it anisotropic later...
    InterpStencil<T> ret(Point::Ones(a_ratio));
    Box K = Box(Point::Ones(a_ratio-1));
    Box U = Box::Cube(3).shift(Point::Ones(-1));

    double r = a_ratio;
     
    for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
    {
        Point n = *kiter;

        Stencil<double>& S = ret(*kiter);
        for (auto iter = U.begin(); iter != U.end(); ++iter)
        {
            T coef = 1;
            for (int ii = 0; ii < DIM; ii++)
            {
                coef *= (T)qcoef(n[ii],r,(*iter)[ii]);
            }
            S += coef*Shift(*iter);
        }
        S.destRatio() = Point::Ones(a_ratio);
        S.destShift() = n;
    }
    return ret;    
}

// helper function for Build(...)
/*
    Returns <x^p> over the coarse cell shifted by s
*/
inline double coarseCoef(Point p, Point s, double h)
{
    double coef = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        double high = h*(2.0*s[d] + 1)/2.0;
        double low = h*(2.0*s[d] - 1)/2.0;
        coef *= 1.0/(p[d]+1)
            *(std::pow(high, p[d]+1) - std::pow(low, p[d]+1));
    }
    coef /= std::pow(h,DIM);
    return coef;
}

// helper function for Build(...)
/**
    Returns <x^p> over the fine cell s with coarse spacing h
*/
inline double fineCoef(Point p, Point s, double h, int refRatio)
{
    double hf = h/refRatio;
    auto low = Proto::Point::Ones(-refRatio/2) + s;
    auto high = low + Proto::Point::Ones();
    
    double coef = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        coef *= 1.0/(p[d]+1)
            *(std::pow(high[d]*hf, p[d]+1) - std::pow(low[d]*hf, p[d]+1));
    }
    coef /= std::pow(hf,DIM);
    return coef;
}

// LAPACK functions for Build
extern "C" {
    int dgemv_(char*, int*, int*, double*, double*,
            int*, double*, int*, double*, double*, int*);
    int dgemm_(char*, char*, int*, int*, int*, double*,
            double*, int*, double*, int*, double*, double*, int*);
    int dgetri_(int*, double*, int*, int*, double*, int*, int*);
    int dgetrf_(int*, int*, double*, int*, int*, int*);
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(int a_smax, Box a_skernel,
        int a_expMax, int a_refRatio)
{
    std::vector<Point> shifts;
    for (auto iter = a_skernel.begin(); iter != a_skernel.end(); ++iter)
    {
        int dist = (*iter).absSum();
        if (dist <= a_smax && dist > 0)
        {
            shifts.push_back(*iter);
        }
    }
    std::vector<Point> orders;
    Box pkernel = Box::Cube(a_expMax+1);
    for (auto iter = pkernel.begin(); iter != pkernel.end(); ++iter)
    {
        int sum = (*iter).sum();
        if (sum <= a_expMax){orders.push_back(*iter);}
    }
    return InterpStencil<T>::Build(shifts,orders,a_refRatio);
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(
        std::vector<Point>& a_shifts, 
        int a_maxOrder,
        int a_refRatio)
{
    std::vector<Point> orders;
    Box K = Box(Point::Ones(a_maxOrder));
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        if ((*iter).sum() <= a_maxOrder)
        {
            orders.push_back(*iter);
        }
    }

    return InterpStencil<T>::Build(a_shifts, orders, a_refRatio);
}

template<class T>
InterpStencil<T> InterpStencil<T>::Build(
        std::vector<Point>& a_shifts, 
        const std::vector<Point>& a_orders,
        int a_refRatio)
{
    // ==============================================================
    // COMPUTE COEFFICIENTS
    // --------------------
    
    // Prevent the zero shift or duplicate shifts from being used
  using std::set;
    std::set<Point> tmp;
    for (auto iter = a_shifts.begin(); iter != a_shifts.end(); ++iter)
    {
        if ((*iter) != Proto::Point::Zeros())
        {
            tmp.insert(*iter);
        }
    }
    a_shifts.clear();
    for (auto iter = tmp.begin(); iter != tmp.end(); ++iter)
    {
        a_shifts.push_back(*iter);
    }
    int M = a_shifts.size();
    int N = a_orders.size();
    /*
    std::cout << "Shifts: " << std::endl;
    for (auto iter = a_shifts.begin(); iter != a_shifts.end(); ++iter)
    {
        std::cout << "\t" << *iter << std::endl;
    }
    std::cout << "Orders: " << std::endl;
    for (auto iter = a_orders.begin(); iter != a_orders.end(); ++iter)
    {
        std::cout << "\t" << *iter << std::endl;
    }
    */
    double* X = (double*)malloc(sizeof(double)*M*N);
    for (int ii = 0; ii < M; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            X[jj*M+ii] = coarseCoef(a_orders[jj], a_shifts[ii], 1.0);
        }
    }
    
    double* X0 = (double*)malloc(sizeof(double)*N);
    for (int jj = 0; jj < N; jj++)
    {
        X0[jj] = coarseCoef(a_orders[jj],Point::Zeros(),1.0);
    }
    
    double* XTX2 = (double*)malloc(sizeof(double)*N*N);
    char _trans_a, _trans_b;
    int _m, _n, _k, _lda, _ldb, _ldc;
    double _alpha, _beta;

    _trans_a = 'T'; _trans_b = 'N';
    _m = N; _n = N; _k = M;
    _lda = M; _ldb = M; _ldc = N;
    _alpha = 2.0; _beta = 0.0;
    dgemm_(&_trans_a, &_trans_b, &_m, &_n, &_k, &_alpha,
            X, &_lda, X, &_ldb, &_beta, XTX2, &_ldc);
    
    double* Y = (double*)malloc(sizeof(double)*(N+1)*(N+1));
    for (int ii = 0; ii <= N; ii++)
    {
        for (int jj = 0; jj <= N; jj++)
        {
            if (ii < N && jj < N) {
                Y[jj*(N+1)+ii] = XTX2[jj*N+ii];
            } else if (ii >= N && jj < N){
                Y[jj*(N+1)+ii] = X0[jj];
            } else if (jj >= N && ii < N){
                Y[jj*(N+1)+ii] = X0[ii];
            } else {
                Y[jj*(N+1)+ii] = 0;
            }
        }
    }

    double* Z = (double*)malloc(sizeof(double)*(N+1)*(M+1));
    for (int ii = 0; ii <= N; ii++)
    {
        for (int jj = 0; jj <= M; jj++)
        {
            if (ii < N && jj < M) {
                Z[jj*(N+1)+ii] = 2.0*X[ii*M+jj];//X transpose
            } else if (ii >= N && jj < M){
                Z[jj*(N+1)+ii] = 0;
            } else if (jj >= M && ii < N){
                Z[jj*(N+1)+ii] = 0;
            } else {
                Z[jj*(N+1)+ii] = 1;
            }
        }
    }

    int* pivots = (int*)malloc(sizeof(int)*(N+1));
    double* _work = (double*)malloc(sizeof(double)*(N+1));
    int _info;
    _n = N+1,
    dgetrf_(&_n, &_n, Y, &_n, pivots, &_info);
    PROTO_ASSERT(_info == 0, "Y is singular.");
    dgetri_(&_n, Y, &_n, pivots, _work, &_n, &_info);
    PROTO_ASSERT(_info == 0, "Failed to invert Y.");
    
    double* YinvZ = (double*)malloc(sizeof(double)*(N+1)*(M+1));
    _trans_a = 'N'; _trans_b = 'N';
    _n = M+1; _k = N+1; _m = N+1;
    _alpha = 1.0;
    _lda = _m; _ldb = _m; _ldc = _m;
    dgemm_(&_trans_a, &_trans_b, &_m, &_n, &_k,
            &_alpha, Y, &_lda, Z, &_ldb, &_beta, YinvZ, &_ldc);

    double* Q = (double*)malloc(sizeof(double)*N*(M+1));
    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < M+1; jj++)
        {
            Q[N*jj+ii] = YinvZ[(N+1)*jj+ii];
        }
    }
   
    int Mf = ipow<DIM>(a_refRatio);
    auto fkernel = Box::Cube(a_refRatio);
    double* Xf = (double*)malloc(sizeof(double)*Mf*N);
    int ii = 0;
    for (auto fiter = fkernel.begin(); fiter != fkernel.end();
             ++fiter, ++ii)
    {
        for (int jj = 0; jj < N; jj++)
        {
            Xf[jj*Mf+ii] = fineCoef(a_orders[jj], *fiter, 1.0, a_refRatio);
        }
    }
    
    double* R = (double*)malloc(sizeof(double)*Mf*(M+1));
    _trans_a = 'N'; _trans_b = 'N';
    _n = M+1; 
    _k = N; _m = Mf;
    _alpha = 1.0; _beta = 0.0;
    _lda = _m; _ldb = _k; _ldc = _m;
    dgemm_(&_trans_a, &_trans_b, &_m, &_n, &_k, &_alpha,
            Xf, &_lda, Q, &_ldb, &_beta, R, &_ldc);

    free(X);
    free(X0);
    free(XTX2);
    free(Y);
    free(Z);
    free(YinvZ);
    free(_work);
    free(pivots);
    free(Q);
    free(Xf);
    
    // ==============================================================
    // BUILD STENCILS
    // --------------

    InterpStencil<T> ret(Point::Ones(a_refRatio));
    ii = 0;
    for (auto fiter = fkernel.begin(); fiter != fkernel.end();
             ++fiter, ++ii)
    {
        Stencil<T>& S = ret(*fiter);
        for (int jj = 0; jj < a_shifts.size(); jj++)
        {
            S += ((T)R[Mf*jj+ii])*Shift(a_shifts[jj]);
        }
        S += ((T)R[Mf*M+ii])*Shift::Zeros();
        S.destRatio() = Point::Ones(a_refRatio);
        S.destShift() = (*fiter);
    }
    return ret;
}

#endif
