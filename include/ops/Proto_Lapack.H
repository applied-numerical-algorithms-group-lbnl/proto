#pragma once
#ifndef _PROTO_LAPACK_H_
#define _PROTO_LAPACK_H_

#define PROTO_LAPACK(NAME, name) d ## name ## _

/// Forward declaration of lapack routines
extern "C" 
{
    /// Double General Solve
    void PROTO_LAPACK(GESV,gesv)(int*, int*, double*, int*, int*, double*, int*, int*);

    /// Double Least Squares 
    void PROTO_LAPACK(GELS,gels)(char* TRANS, int* M, int* N, int* NRHS, double* A, 
            int* LDA, double* B, int* LDB, 
            double* WORK, int* LWORK, int* INFO);

    /// Double Least Squares transpose via SVD
    void PROTO_LAPACK(GELSD,gelsd)(int* N, int* M, int* NRHS, double* At, 
            int* LDA, double* B, int* LDB, double* S,
            double* RCOND, int* RANK, 
            double* WORK, int* LWORK, int* IWORK, int* INFO);

#ifndef PR_BLIS
    /// Matrix multiply
    void PROTO_LAPACK(GEMV,gemv)(char* TRANS, int* M, int* N, double* ALPHA, 
            double* A, int* LDA, double* X, int* INCX,
            double* BETA, double* Y, int* INCY);
#endif
    /// Double Equality-Constrained Least Squares 
    void PROTO_LAPACK(GGLSE,gglse)(int* M, int* N, int* P, double* A, int* LDA, double* B,
            int* LDB, double* C, double* D, double* X, double* WORK, int* LWORK, int* INFO);

    /// QR factorization with Pivoting
    void PROTO_LAPACK(GEQP3,geqp3)(int* M, int* N, double* A, int* LDA, int* JPVT,
            double* TAU, double* WORK, int* LWORK, int* INFO);

    /// QR factorization
    void PROTO_LAPACK(GEQRF,geqrf)(int* M, int* N, double* A, int* LDA, 
            double* TAU, double* WORK, int* LWORK, int* INFO);

    /// Multiply by Q from QR factorization
    void PROTO_LAPACK(ORMQR,ormqr)(char* SIDE, char* TRANS, int* M, int* N, int* K, 
            double* A, int* LDA, double* TAU, double* C, int* LDC,
            double* WORK, int* LWORK, int* INFO);

    /// Generate Q from QR factorization
    void PROTO_LAPACK(ORGQR,orgqr)(int* M, int* N, int* K, double* A, int* LDA, 
            double* TAU, double* WORK, int* LWORK, int* INFO);

    /// Triangular solve
    void PROTO_LAPACK(TRTRS,trtrs)(char* UPLO, char* TRANS, char* DIAG, int* N, int* NRHS, 
            double* A, int* LDA, double* B, int* LDB, int* INFO);

    /// LU decomoposition of a general matrix
    void PROTO_LAPACK(GETRF,getrf)(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    /// generate inverse of a matrix given its LU decomposition
    void PROTO_LAPACK(GETRI,getri)(int* N, double* A, int* lda,
            int* IPIV, double* WORK, int* lwork, int* INFO);

    /// get condition number of the matrix
    void PROTO_LAPACK(GECON,gecon)(char* TRANS, int* N,  double* A, int* LDA, double* ANORM,
            double* rcond, double* WORK, int* LWORK, int* INFO);

    /// get condition number of triangular matrix
    void PROTO_LAPACK(TRCON,trcon)(const char* TRANS, const char* UPLO, const char* DIAG,
            const int* N,  const double* A, const int* LDA,
            double* rcond, double* WORK, int* IWORK, int* INFO);


    /// multiply two matrices (BLAS routine)
    void PROTO_LAPACK(GEMM,gemm)(const char* TRANSA, const char* TRANSB,
            const int* M, const int* N, const int* K, const double* ALPHA,
            const double* A, const int* LDA,
            const double* B, const int* LDB,
            const double* BETA, double* C, const int* LDC);

    /// solve matrix equations with triangular matrix (BLAS routine)
    void PROTO_LAPACK(TRSM,trsm)(const char* SIDE, const char* UPLO,
            const char* TRANSA, const char* DIAG,
            const int* M, const int* N, const double* ALPHA,
            const double* A, const int* LDA,
            double* B, const int* LDB);
}
#endif //end include guard
