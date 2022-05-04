CoordPermutation::CoordPermutation()
{
    for (int ii = 0; ii < DIM; ii++)
    {
        m_matrix[ii].fill(0);
        m_matrix[ii][ii] = 1;
    }
}

CoordPermutation::CoordPermutation(
        std::initializer_list<std::array<int, 3>> a_mapping)
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
    for (auto map : a_mapping)
    {
        m_matrix[map[1]].fill(0);
        m_matrix[map[1]][map[0]] = map[2];
    }
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

CoordPermutation CoordPermutation::cw()
{
#if DIM >= 2
    return CoordPermutation{{0,1,1},{1,0,-1}};
#else
    return CoordPerutation();
#endif
}

CoordPermutation CoordPermutation::ccw()
{
#if DIM >= 2
    return CoordPermutation{{0,1,-1},{1,0,1}};
#else
    return CoordPerutation();
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

void CoordPermutation::print() const
{
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
    for (auto row : m_matrix)
    {
        for (auto val : row)
        {
            pout() << val << ", ";
        }
        pout() << std::endl;
    }
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

