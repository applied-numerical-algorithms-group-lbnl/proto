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
            EXPECT_EQ(M(ii,jj), 17);
        }
    }
    initialize(M);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            double val = jj*M.M() + ii;
            EXPECT_EQ(M(ii,jj), val);
        }
    }
    Matrix<double> A{
        {0, 4, 8},
        {1, 5, 9},
        {2, 6, 10},
        {3, 7, 11}};
    EXPECT_EQ(M, A);
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
            EXPECT_EQ(M(ii,jj), val);
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
            EXPECT_EQ(M(ii,jj), val);
        }
    }
}
#ifdef PR_BLIS
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
            EXPECT_EQ(Ai(0, jj), val);
        }
        Ai.set(ii);
        for (int jj = 0; jj < n; jj++)
        {
            EXPECT_EQ(A(ii, jj), ii);
        }

    }
    for (int jj = 0; jj < n; jj++)
    {
        auto Aj = A.col(jj);
        for (int ii = 0; ii < m; ii++)
        {
            EXPECT_EQ(Aj(ii, 0), ii);
        }
        Aj.set(jj);
        for (int ii = 0; ii < m; ii++)
        {
            EXPECT_EQ(A(ii, jj), jj);
        }
    }
}
#endif
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
            EXPECT_EQ(A(kk), A(ii, jj));
        }
    }
    for (int kk = 0; kk < m*n; kk++)
    {
        double val = m*n - kk;
        int ii = kk % m;
        int jj = kk / m;
        A.set(kk, val);
        EXPECT_EQ(A(ii, jj), val);
    }
}

TEST(Matrix, Copy)
{
    int m = 4;
    int n = 3;
    Matrix<double> A(m, n);
    Matrix<double> B(m, n);
    initialize(A);
    B.set(0);
    A.copyTo(B);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            EXPECT_EQ(A(ii,jj), B(ii,jj));
        }
    }
}

TEST(Matrix, TransposeCopy)
{
    int m = 4;
    int n = 3;
    Matrix<double> A(m, n);
    Matrix<double> B(n, m);
    initialize(A);
    B.set(0);
    A.transpose().copyTo(B);
    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            EXPECT_EQ(A(ii,jj), B(jj,ii));
        }
    }
}

TEST(Matrix, AddSubtract)
{
    int m = 4;
    int n = 3;
    Matrix<double> A(m, n);
    Matrix<double> B(m, n);
    Matrix<double> C(m, n);
    Matrix<double> D(m, n);
    Matrix<double> I(m, n);
    Matrix<double> J(m, n);
    initialize(A);
    initialize(B);
    C.set(17);
    D.set(17);
    C += A;
    D -= A;
    auto E = A + B;
    auto F = A - B;
    auto G = A + 1;
    auto H = A - 1;
    initialize(I);
    initialize(J);
    I += 1;
    J -= 1;

    for (int ii = 0; ii < m; ii++)
    {
        for (int jj = 0; jj < n; jj++)
        {
            int kk = jj*m+ii;
            EXPECT_EQ(A(ii,jj), kk);
            EXPECT_EQ(B(ii,jj), kk);
            EXPECT_EQ(C(ii,jj), 17+kk);
            EXPECT_EQ(D(ii,jj), 17-kk);
            EXPECT_EQ(E(ii,jj), 2*kk);
            EXPECT_EQ(F(ii,jj), 0);
            EXPECT_EQ(G(ii,jj), kk+1);
            EXPECT_EQ(H(ii,jj), kk-1);
            EXPECT_EQ(I(ii,jj), kk+1);
            EXPECT_EQ(J(ii,jj), kk-1);
        }
    }
}

TEST(Matrix, Multiply)
{
    int m = 4;
    int n = 3;
    int p = 3;
    int q = 4;
    Matrix<double> A(m, n);
    Matrix<double> B(p, q);
    Matrix<double> C(m, n);
    initialize(A);
    initialize(B);
    initialize(C);
    auto AB = A*B;
    auto BA = B*A;
    
    C *= 17;
    auto D = A*17;

    EXPECT_EQ(AB.M(), m);
    EXPECT_EQ(AB.N(), q);
    EXPECT_EQ(BA.M(), p);
    EXPECT_EQ(BA.N(), n);

    for (int ii = 0; ii < m; ii++)
    for (int jj = 0; jj < n; jj++)
    {
        EXPECT_EQ(C(ii,jj), A(ii,jj)*17);
        EXPECT_EQ(D(ii,jj), A(ii,jj)*17);
    }
    for (int ii = 0; ii < m; ii++)
    for (int jj = 0; jj < q; jj++)
    {
        double Cij = 0;
        for (int kk = 0; kk < n; kk++)
        {
            Cij += A(ii,kk)*B(kk,jj);
        }
        EXPECT_EQ(AB(ii,jj), Cij);
    }
    for (int ii = 0; ii < p; ii++)
    for (int jj = 0; jj < n; jj++)
    {
        double Cij = 0;
        for (int kk = 0; kk < m; kk++)
        {
            Cij += B(ii,kk)*A(kk,jj);
        }
        EXPECT_EQ(BA(ii,jj), Cij);
    }
}

TEST(Matrix, Transpose)
{
    int m = 4;
    int n = 3;
    Matrix<double> A(m,n);
    initialize(A);
    auto AT = A.transpose();
    auto ATA = A.transpose()*A;

    Matrix<double> ATcopy(n,m);
    A.transpose().copyTo(ATcopy);

    EXPECT_EQ(AT.M(), n);
    EXPECT_EQ(AT.N(), m);
    EXPECT_EQ(ATcopy.M(), n);
    EXPECT_EQ(ATcopy.N(), m);
    EXPECT_EQ(ATA.M(), n);
    EXPECT_EQ(ATA.N(), n);
    for (int ii = 0; ii < m; ii++)
    for (int jj = 0; jj < n; jj++)
    {
        double Aij = m*jj + ii;
        EXPECT_EQ(A(ii,jj), Aij);
        EXPECT_EQ(A.get(ii,jj), Aij);
        EXPECT_EQ(AT(jj, ii), Aij);
        EXPECT_EQ(AT.get(jj, ii), Aij);
        EXPECT_EQ(ATcopy(jj, ii), Aij);
    }
    for (int ii = 0; ii < n; ii++)
    for (int jj = 0; jj < n; jj++)
    {
        double val = 0;
        for (int kk = 0; kk < m; kk++)
        {
            val += AT(ii,kk)*A(kk,jj);
        }
        EXPECT_EQ(ATA(ii,jj), val);
    }
}

TEST(Matrix, Inverse)
{
    int m = 2;
    Matrix<double> A(m,m);
    initialize(A);
    A += 1;
    auto B = A.inverse();
    auto C = A*B;
    auto I = Matrix<double>::I(m);
    for (int kk = 0; kk < m*m; kk++)
    {
        EXPECT_NEAR(C(kk), I(kk), 1e-12);
    }
}

TEST(Matrix, LeastSquares)
{
    Matrix<double> A{
        {2, -1},
        {1,  2},
        {1,  1}};
    Matrix<double> b{{2},{1},{4}};
    Matrix<double> x{{10.0/7.0},{3.0/7.0}};
    
    solveLS(A, b);

    double error = pow(b(0,0) - x(0,0), 2) + pow(b(1,0) - x(1,0), 2);
    EXPECT_LT(error, 1e-16);
}

TEST(Matrix, SVD)
{
    Matrix<double> A{
        { 2, -1,  0},
        { 4,  3, -2}};
    
    Matrix<double> U, S, VT;
    
    A.svd(U, S, VT);
   
    auto Sigma = Matrix<double>::I(S.M(), S.M());
    for (int ii = 0; ii < S.M(); ii++)
    {
        Sigma(ii,ii) = S(ii,0);
    }

    auto B = U*(Sigma*VT);
    B -= A;
    double error = B.normInf();
    EXPECT_LT(B.normInf(), 1e-12);
}

TEST(Matrix, PseudoInverse)
{
    Matrix<double> A{
        {1, -1},
        {0,  1},
        {1,  0}};

    Matrix<double> Soln{
        { 1, 1, 2},
        {-1, 2, 1}};
    
    Soln *= (1.0/3);

    auto B = A.pseudoInverse();

    auto ABA = A*(B*A);
    auto BAB = B*(A*B);

    ABA -= A;
    BAB -= B;
    Soln -= B;
    EXPECT_LT(Soln.normInf(), 1e-12);
    EXPECT_LT(ABA.normInf(), 1e-12);
    EXPECT_LT(BAB.normInf(), 1e-12);
}
int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
#ifdef PR_BLIS
    cout << "RUNNING MATRIX TESTS WITH BLIS BACKEND" << std::endl; 
#else
    cout << "RUNNING MATRIX TESTS WITHOUT BLIS BACKEND" << std::endl; 
#endif
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
