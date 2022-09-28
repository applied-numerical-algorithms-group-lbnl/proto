#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

template<typename T>
void initialize(Matrix<T>& a_A)
{
    for (int ii = 0; ii < a_A.M(); ii++)
    {
        for (int jj = 0; jj < a_A.N(); jj++)
        {
            T val = jj*a_A.M() + ii;
            a_A.set(ii, jj, val);
        }
    }
}

TEST(Matrix, Construction) {
    int m = 4;
    int n = 3;
    Matrix<double> M(m, n);
    EXPECT_EQ(M.numRows(), m);
    EXPECT_EQ(M.numCols(), n);
    M.set(17);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            EXPECT_EQ(M.get(ii,jj), 17);
        }
    }
    initialize(M);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            double val = jj*M.M() + ii;
            EXPECT_EQ(M.get(ii,jj), val);
        }
    }
}

TEST(Matrix, BufferConstruction)
{
    int m = 4;
    int n = 3;
    std::shared_ptr<double> buffer((double*)proto_malloc<HOST>(m*n*sizeof(double)));
    for (int ii = 0; ii < m*n; ii++)
    {
        *(buffer.get()+ii) = ii;
    }
    Matrix<double> M(buffer, m, n);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            double val = jj*M.M() + ii;
            EXPECT_EQ(M.get(ii,jj), val);
        }
    }
    for (int ii = 0; ii < m*n; ii++)
    {
        *(buffer.get()+ii) = (m*n - ii);
    }
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            double val = jj*M.M() + ii;
            val = m*n - val;
            EXPECT_EQ(M.get(ii,jj), val);
        }
    }
}

TEST(Matrix, Slicing)
{
    int m = 4;
    int n = 3;
    Matrix<double> A(m, n);
    initialize(A);
    for (int ii = 0; ii < m; ii++)
    {
        auto Ai = A.row(ii);
        for (int jj = 0; jj < n; jj++)
        {
            double val = jj*m + ii;
            EXPECT_EQ(Ai.get(0, jj), val);
        }
        Ai.set(ii);
        for (int jj = 0; jj < n; jj++)
        {
            EXPECT_EQ(A.get(ii, jj), ii);
        }

    }
    for (int jj = 0; jj < n; jj++)
    {
        auto Aj = A.col(jj);
        for (int ii = 0; ii < m; ii++)
        {
            EXPECT_EQ(Aj.get(ii, 0), ii);
        }
        Aj.set(jj);
        for (int ii = 0; ii < m; ii++)
        {
            EXPECT_EQ(A.get(ii, jj), jj);
        }
    }
}

TEST(Matrix, LinearIndexing)
{
    int m = 4;
    int n = 3;
    Matrix<double> A(m, n);
    initialize(A);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            int kk = jj*m+ii;
            EXPECT_EQ(A.get(kk), A.get(ii, jj));
        }
    }
    for (int kk = 0; kk < m*n; kk++)
    {
        double val = m*n - kk;
        int ii = kk % m;
        int jj = kk / m;
        A.set(kk, val);
        EXPECT_EQ(A.get(ii, jj), val);
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
