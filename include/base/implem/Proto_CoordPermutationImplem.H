CoordPermutation::CoordPermutation()
{
    for (int ii = 0; ii < DIM; ii++)
    {
        m_matrix[ii].fill(0);
        m_matrix[ii][ii] = 1;
    }
}
/*
CoordPermutation::CoordPermutation(std::array<std::array<int, DIM>, DIM> a_matrix)
{
    m_matrix = a_matrix;
}
*/

CoordPermutation::CoordPermutation(
        std::vector<std::array<int, 3>>& a_mapping)
{
    define(a_mapping);
}

CoordPermutation::CoordPermutation(
        std::initializer_list<std::array<int, 3>> a_mapping)
{
    std::vector<std::array<int, 3>> v(a_mapping);
    define(v);
}

CoordPermutation::CoordPermutation(
        std::initializer_list<std::array<Point, 2>> a_mapping)
{
    std::vector<std::array<int, 3>> mapping;
    for (auto map : a_mapping)
    {
        int c0, c1, s;
        s = 1;
        for (int dir = 0; dir < DIM; dir++)
        {
            if (map[0][dir] != 0) {c0 = dir; s *= map[0][dir]; }
            if (map[1][dir] != 0) {c1 = dir; s *= map[1][dir]; }
        }
        std::array<int, 3> arr{c0, c1, s/abs(s)};
        mapping.push_back(arr);
    }
    define(mapping);
}

void CoordPermutation::define(std::vector<std::array<int, 3>>& a_mapping)
{
    for (auto map : a_mapping)
    {
        PROTO_ASSERT(map[0] >= 0 && map[0] < DIM,
                "CoordPermutation::Constructor | Invalid mapping input: [%i, %i, %i].",
                map[0], map[1], map[2]);
        PROTO_ASSERT(map[1] >= 0 && map[1] < DIM,
                "CoordPermutation::Constructor | Invalid mapping input: [%i, %i, %i].",
                map[0], map[1], map[2]);
        PROTO_ASSERT(map[2] == 1 || map[2] == -1,
                "CoordPermutation::Constructor | Invalid mapping input: [%i, %i, %i].",
                map[0], map[1], map[2]);
    }
    for (int ii = 0; ii < DIM; ii++)
    {
        m_matrix[ii].fill(0);
        m_matrix[ii][ii] = 1;
    }
    for (auto map : a_mapping)
    {
        m_matrix[map[1]].fill(0);
        m_matrix[map[1]][map[0]] = map[2];
    }
    PROTO_ASSERT(this->isValid(),
            "CoordPermutation::Constructor | Resulting object is not a valid permutation.");
}

void CoordPermutation::defineMatrix(std::array<std::array<int, DIM>, DIM> a_matrix)
{
    m_matrix = a_matrix;
    PROTO_ASSERT(this->isValid(),
            "CoordPermutation::Constructor | Resulting object is not a valid permutation.");
}

CoordPermutation CoordPermutation::identity()
{
    return CoordPermutation();
}

CoordPermutation CoordPermutation::reverse()
{
    CoordPermutation ret;
    for (int ii = 0; ii < DIM; ii++)
    {
        ret.m_matrix[ii][ii] = -1;
    }
    return ret;
}

CoordPermutation CoordPermutation::cw(int a_coord)
{
#if DIM == 2
    //return cw(2);
    return CoordPermutation{{0,1,1},{1,0,-1}};
#else
    int c0 = ((a_coord - 2) + 3) % 3;
    int c1 = ((a_coord - 1) + 3) % 3;
    return CoordPermutation{{c0, c1,1},{c1,c0,-1}};
#endif
}

CoordPermutation CoordPermutation::rotatePlane(int a_c0, int a_c1)
{
#if DIM >= 2
    return CoordPermutation({{a_c0,a_c1,1},{a_c1,a_c0,-1}});
#else
    return CoordPermutation();
#endif
}

CoordPermutation CoordPermutation::ccw(int a_coord)
{
#if DIM == 2
    //return ccw(2);
    return CoordPermutation{{0,1,-1},{1,0,1}};
#else
    int c0 = ((a_coord - 2) + 3) % 3;
    int c1 = ((a_coord - 1) + 3) % 3;
    return CoordPermutation{{c0, c1,-1},{c1,c0,1}};
#endif
}

Point CoordPermutation::operator()(Point a_pt) const
{
    Point mapped;
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (m_matrix[jj][ii] != 0)
            {
                mapped[jj] = m_matrix[jj][ii]*a_pt[ii];
            }
        }
    }
    return mapped;
}

CoordPermutation CoordPermutation::inverse() const
{
    std::array<std::array<int, DIM>, DIM> invMatrix;
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            invMatrix[ii][jj] = m_matrix[jj][ii];
        }
    }
    CoordPermutation R;
    R.defineMatrix(invMatrix);
    return R;
}

Point CoordPermutation::rotatePoint(Point a_point, const Box& a_srcBox, const Box& a_dstBox) const
{
    //PROTO_ASSERT(a_box.low() == Point::Zeros(),
    //        "CoordPermutation::rotatePoint | Error: input Box must have an origin at (0,...,0)");
    Point src0 = a_srcBox.low();
    Point dst0 = a_dstBox.low();
    Box B0 = a_srcBox.shift(-src0);
    Point a = -(B0.sizes());
    Point b = (*this)(a);
    Point s;
    for (int d = 0; d < DIM; d++)
    {
        s[d] = max(0,b[d]);
    }
    Point r = (*this)(a_point - src0) + s + dst0;
    //std::cout << "Rotated " << a_point << " in box " << a_box << " and got " << r << std::endl;
    return r;
}

Point CoordPermutation::rotateCell(Point a_point, const Box& a_srcBox, const Box& a_dstBox) const
{
    //PROTO_ASSERT(a_box.low() == Point::Zeros(),
    //        "CoordPermutation::rotatePoint | Error: input Box must have an origin at (0,...,0)");
    Point src0 = a_srcBox.low();
    Point dst0 = a_dstBox.low();
    Box B0 = a_srcBox.shift(-src0);
    Point a = -(B0.sizes() - Point::Ones());
    Point b = (*this)(a);
    Point s;
    for (int d = 0; d < DIM; d++)
    {
        s[d] = max(0,b[d]);
    }
    Point r = (*this)(a_point - src0) + s + dst0;
    //std::cout << "Rotated " << a_point << " in box " << a_box << " and got " << r << std::endl;
    return r;
}

template< typename T, unsigned int C, unsigned char D, unsigned char E>
void CoordPermutation::rotateBuffer(
        T* srcData, T* dstData,
        const Box& srcBox, const Box& dstBox) const
{
    MemType srcMem = pointerMemType(srcData);
    MemType dstMem = pointerMemType(srcData);
    PROTO_ASSERT(srcMem == HOST,
            "rotateBuffer | Error: source data buffer is not HOST.");
    PROTO_ASSERT(dstMem == HOST,
            "rotateBuffer | Error: destination data buffer is not HOST.");
    PROTO_ASSERT(srcBox.size() == dstBox.size(),
            "rotateBuffer | Error: \
            Rotated box must be the same size as the current box.");
    unsigned int N = srcBox.size();
    unsigned int L = C*D*E*N;
    for (unsigned int kk = 0; kk < srcBox.size(); kk++)
    {
        unsigned int ii = kk;
        Point pi = srcBox[ii];
        Point pj = this->rotateCell(pi, srcBox, dstBox);
        unsigned int jj = dstBox.index(pj);
        
        for (int ee = 0; ee < E; ee++)
        for (int dd = 0; dd < D; dd++)
        for (int cc = 0; cc < C; cc++)
        {
            PROTO_ASSERT(ii < L,
                    "BoxData::rotate | Error: unrotated index out of bounds \
                    (%u > %u)", ii, L);
            PROTO_ASSERT(jj < L,
                    "BoxData::rotate | Error: rotated index out of bounds \
                    (%u > %u)", jj, L);
            dstData[jj] = srcData[ii];
            ii += N;
            jj += N;
        }
    }
}

bool CoordPermutation::operator==(const CoordPermutation& a_rhs) const
{
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (m_matrix[ii][jj] != a_rhs.m_matrix[ii][jj]) { return false; }
        }
    }
    return true;
}

void CoordPermutation::print() const
{
    using Proto::pout;
    for (int jj = 0; jj < DIM; jj++)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            if (m_matrix[jj][ii] != 0)
            {
                pout() << ii << " --> " << jj;
                pout() << " (" << m_matrix[jj][ii] << ")" << std::endl;
            }
        }
    }
    this->printMatrix();
}

void CoordPermutation::printMatrix() const
{
    using Proto::pout;
    for (auto row : m_matrix)
    {
        for (auto val : row)
        {
            pout() << val << ", ";
        }
        pout() << std::endl;
    }
}

bool CoordPermutation::isIdentity() const
{
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < DIM; jj++)
        {
            if ((ii != jj) && (m_matrix[ii][jj] != 0)) { return false; }
            if ((ii == jj) && (m_matrix[ii][jj] != 1)) { return false; }
        }
    }
    return true;
}

bool CoordPermutation::isValid() const
{
    for (auto row : m_matrix)
    {
        int nonZeros = 0;
        for (auto val : row)
        {
            if (val != 0)
            {
                nonZeros++;
                if (val < -1 || val > 1)
                {
                    return false;
                }
            }
        }
        if (nonZeros != 1) { return false; }
    }
    return true;
}

CoordPermutation operator*(const CoordPermutation& a_A, const CoordPermutation& a_B)
{
    auto A = a_A.matrix();
    auto B = a_B.matrix();
    std::array<std::array<int, DIM>, DIM> M;
    for (int ii = 0; ii < DIM; ii++)
    {
        for (int jj = 0; jj < DIM; jj++)
        {
            int Mij = 0;
            for (int kk = 0; kk < DIM; kk++)
            {
                Mij += A[kk][ii]*B[jj][kk];
            }
            M[jj][ii] = Mij;
        }
    }
    CoordPermutation R;
    R.defineMatrix(M);
    return R;
}

