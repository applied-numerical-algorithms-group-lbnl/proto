
template<typename T>
Matrix<T> Matrix<T>::Identity(int a_M, int a_N)
{
    PR_TIME("Matrix::Identity");
    Matrix<T> ret(a_M, a_N);
    ret.set(0);
    for (int ii = 0; ii < a_M; ii++)
    for (int jj = 0; jj < a_N; jj++)
    {
        if (ii == jj) {ret(ii,jj) = 1;}
        else {ret(ii,jj) = 0;}
    }
    return ret;
}

template<typename T>
Matrix<T>::Matrix(unsigned int a_numRows, unsigned int a_numCols, MatrixStorageType a_storage)
{
    define(a_numRows, a_numCols, a_storage);
}

template<typename T>
Matrix<T>::Matrix(
        std::shared_ptr<T>& a_buffer,
        unsigned int a_numRows,
        unsigned int a_numCols,
        MatrixStorageType a_storage)
{
    PR_TIME("Matrix(buffer, numRows, numCols)");
    m_storage = a_storage;
#ifdef PR_BLIS
    int row_stride = 1;
    int col_stride = a_numRows;
    if (a_storage == ROW_MAJOR)
    {
        row_stride = a_numCols;
        col_stride = 1;
    }
    bli_obj_create_without_buffer(blisDatatype(), a_numRows, a_numCols, &m_object);
    bli_obj_attach_buffer(a_buffer.get(), row_stride, col_stride, 0, &m_object);
#else
    m_transpose = false;
    m_numRows = a_numRows;
    m_numCols = a_numCols;
    m_data = a_buffer;
#endif
    m_alias = true;
}

template<typename T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> a_args)
{
    PR_TIME("Matrix(initializer_list)");
    int M = a_args.size();
    int N = (*a_args.begin()).size();
    define(M,N);
    int index = 0;
    int ii = 0;
    int jj = 0;
    for (auto row : a_args)
    {
        PROTO_ASSERT(row.size() == N,
                "Matrix(initializer_list) | Error: inconsistent number of columns specified");
        jj = 0;
        for (auto val : row)
        {
            set(ii, jj, val);
            jj++;
        }
        ii++;
    }
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>&& a_matrix)
{
    PR_TIME("Matrix(Matrix&&)");
    m_alias = a_matrix.m_alias;
    m_storage = a_matrix.m_storage;
#if PR_BLIS
    m_object = std::move(a_matrix.m_object);
    a_matrix.m_alias = true; //prevent double free when rvalue is destroyed 
#else
    m_transpose = a_matrix.m_transpose;
    m_data = a_matrix.m_data;
    m_numRows = a_matrix.m_numRows;
    m_numCols = a_matrix.m_numCols;
#endif
}

template<typename T>
Matrix<T>::~Matrix()
{
    PR_TIME("~Matrix()");
    if (!m_alias)
    {
#ifdef PR_BLIS
        bli_obj_free(&m_object);
#endif
    }
}

template<typename T>
void Matrix<T>::define(
        unsigned int a_numRows,
        unsigned int a_numCols,
        MatrixStorageType a_storage)
{
    PR_TIME("Matrix(numRows, numCols)");
    // using column-major order and suppressing BLIS' buffer packing
#ifdef PR_BLIS
    int row_stride = 1;
    int col_stride = a_numRows;
    if (a_storage == ROW_MAJOR)
    {
        row_stride = a_numCols;
        col_stride = 1;
    }
    bli_obj_create(blisDatatype(), a_numRows, a_numCols, row_stride, col_stride, &m_object);
#else
    m_transpose = false;
    m_numRows = a_numRows;
    m_numCols = a_numCols;
    m_data = std::shared_ptr<T>((T*)proto_malloc<HOST>(a_numRows*a_numCols*sizeof(T)));
#endif
    m_alias = false;
    m_storage = a_storage;
}


template<typename T>
Matrix<T>&
Matrix<T>::operator=(const Matrix<T>&& a_matrix)
{
    PR_TIME("Matrix::operator=(Matrix&&)");
    m_alias = a_matrix.m_alias;
    m_storage = a_matrix.m_storage;
#ifdef PR_BLIS
    m_object = std::move(a_matrix.m_object);
    a_matrix.m_alias = true; //to prevent trying to free the object
#else
    m_transpose = a_matrix.m_transpose;
    m_data = a_matrix.m_data;
    m_numRows = a_matrix.m_numRows;
    m_numCols = a_matrix.m_numCols;
#endif
    return *this;
}

template<typename T>
T 
Matrix<T>::get(unsigned int a_row, unsigned int a_col) const
{
    PR_TIME("Matrix::get(row, col)");
    PROTO_ASSERT(a_row < numRows(), 
            "Matrix::get | Error: row index %u is out of bounds (numRows() = %u)",
            a_row, numRows());
    PROTO_ASSERT(a_col < numCols(), 
            "Matrix::get | Error: col index %u is out of bounds (numCols() = %u)",
            a_col, numCols());
#ifdef PR_BLIS
    T real;
    T imag; //unused
    //TODO: make this work properly if T is complex
    switch (bli_obj_conjtrans_status(&m_object))
    {
        case BLIS_TRANSPOSE:
            bli_getijm(a_col, a_row, &m_object, &real, &imag); break;
        case BLIS_NO_TRANSPOSE:
        default:
            bli_getijm(a_row, a_col, &m_object, &real, &imag); break;
    }
    return real;
#else
    int k = linearIndex(a_row, a_col);
    return m_data.get()[k];
#endif
}

template<typename T>
T
Matrix<T>::get(unsigned int a_index) const
{
    PROTO_ASSERT(a_index < numRows()*numCols(),
            "Matrix::get | Error: index %u is out of bounds", a_index);
    unsigned int ii = a_index % numRows();
    unsigned int jj = a_index / numRows();
    return get(ii, jj);
}

template<typename T>
void 
Matrix<T>::set(unsigned int a_row, unsigned int a_col, const T& a_val)
{
    PR_TIME("Matrix::set(row, col, val)");
    PROTO_ASSERT(a_row < numRows(), 
            "Matrix::get | Error: row index %u is out of bounds (numRows() = %u)",
            a_row, numRows());
    PROTO_ASSERT(a_col < numCols(), 
            "Matrix::get | Error: col index %u is out of bounds (numCols() = %u)",
            a_col, numCols());
    //TODO: make this work properly if T is complex
#ifdef PR_BLIS
    switch (bli_obj_conjtrans_status(&m_object))
    {
        case BLIS_TRANSPOSE:
            bli_setijm(a_val, 0, a_col, a_row, &m_object);
        case BLIS_NO_TRANSPOSE:
        default:
            bli_setijm(a_val, 0, a_row, a_col, &m_object);
    }
#else
    int k = linearIndex(a_row, a_col);
    m_data.get()[k] = a_val;
#endif
}

template<typename T>
void 
Matrix<T>::set(unsigned int a_index, const T& a_val)
{
    PROTO_ASSERT(a_index < numRows()*numCols(),
            "Matrix::get | Error: index %u is out of bounds", a_index);
    unsigned int ii = a_index % numRows();
    unsigned int jj = a_index / numRows();
    return set(ii, jj, a_val);
}

template<typename T>
void 
Matrix<T>::set(const T& a_val)
{
    PR_TIME("Matrix::set(val)");
#ifdef PR_BLIS
    if (a_val == 0)
    {
        bli_setm(&BLIS_ZERO, &m_object);
    } else if (a_val == 1)
    {
        bli_setm(&BLIS_ONE, &m_object);
    } else if (a_val == -1)
    {
        bli_setm(&BLIS_MINUS_ONE, &m_object);
    } else {
        obj_t alpha;
        bli_obj_create_1x1(blisDatatype(), &alpha);
        bli_setsc(a_val, 0, &alpha);
        bli_setm(&alpha, &m_object);
    }
#else
    for (int ii = 0; ii < size(); ii++)
    {
        m_data.get()[ii] = a_val;
    }
#endif
}

template<typename T>
bool
Matrix<T>::operator==(const Matrix<T>& a_rhs) const
{
    PR_TIME("Matrix::operator==");
    if (numRows() != a_rhs.numRows()) { return false; }
    if (numCols() != a_rhs.numCols()) { return false; }

    for (int ii = 0; ii < M(); ii++)
    {
        for (int jj = 0; jj < N(); jj++)
        {
            if (get(ii,jj) != a_rhs.get(ii,jj)) { return false; }
        }
    }
    return true;
}

template<typename T>
bool
Matrix<T>::operator!=(const Matrix<T>& a_rhs) const
{
    return !((*this) == a_rhs);
}

template<typename T>
unsigned int
Matrix<T>::numRows() const
{
    PR_TIME("Matrix::numRows");
#ifdef PR_BLIS
    return bli_obj_length_after_trans(&m_object);
#else
    if (m_transpose) { return m_numCols; }
    return m_numRows;
#endif
}

template<typename T>
unsigned int
Matrix<T>::numCols() const
{
    PR_TIME("Matrix::numCols");
#ifdef PR_BLIS
    return bli_obj_width_after_trans(&m_object);
#else
    if (m_transpose) { return m_numRows; }
    return m_numCols;
#endif
}

#ifndef PR_BLIS
template<typename T>
unsigned int
Matrix<T>::linearIndex(unsigned int a_row, unsigned int a_col) const
{
    switch (m_storage)
    {
        case ROW_MAJOR:
            if (m_transpose)
            {
                return m_numCols*a_col + a_row;
            } else {
                return m_numCols*a_row + a_col;
            }
        case COL_MAJOR:
        default:
            if (m_transpose)
            {
                return m_numRows*a_row + a_col;
            } else {
                return m_numRows*a_col + a_row;
            }
    }
}
#endif

#ifdef PR_BLIS
template<typename T>
Matrix<T>
Matrix<T>::slice(
        unsigned int a_numRows, unsigned int a_numCols,
        unsigned int a_rowOffset, unsigned int a_colOffset)
{
    PR_TIME("Matrix::slice");
    PROTO_ASSERT(a_rowOffset + a_numRows <= numRows(),
            "Matrix::slice | Error: row index is out of bounds");
    PROTO_ASSERT(a_colOffset + a_numCols <= numCols(),
            "Matrix::slice | Error: column index is out of bounds");
    Matrix<T> ret;
    bli_acquire_mpart(
            a_rowOffset, a_colOffset,
            a_numRows, a_numCols,
            &m_object, &(ret.m_object));
    ret.m_alias = true;
    return ret;
}

template<typename T>
const Matrix<T>
Matrix<T>::slice(
        unsigned int a_numRows, unsigned int a_numCols,
        unsigned int a_rowOffset, unsigned int a_colOffset) const
{
    PR_TIME("Matrix::slice");
    PROTO_ASSERT(a_rowOffset + a_numRows <= numRows(),
            "Matrix::slice | Error: row index is out of bounds");
    PROTO_ASSERT(a_colOffset + a_numCols <= numCols(),
            "Matrix::slice | Error: column index is out of bounds");
    Matrix<T> ret;
    bli_acquire_mpart(
            a_rowOffset, a_colOffset,
            a_numRows, a_numCols,
            &m_object, &(ret.m_object));
    ret.m_alias = true;
    return ret;
}

template<typename T>
Matrix<T>
Matrix<T>::row(unsigned int a_rowIndex)
{
    PROTO_ASSERT(a_rowIndex < numRows(),
            "Matrix::row | Error: row index %u is out of bounds (numRows() = %u",
            a_rowIndex, numRows());
    return slice(1, numCols(), a_rowIndex, 0);
}

template<typename T>
Matrix<T>
Matrix<T>::col(unsigned int a_colIndex)
{
    PROTO_ASSERT(a_colIndex < numCols(),
            "Matrix::col | Error: col index %u is out of bounds (numCols() = %u",
            a_colIndex, numCols());
    return slice(numRows(), 1, 0, a_colIndex);
}
#endif

template<typename T>
void
Matrix<T>::copyTo(Matrix<T>& a_rhs) const
{
    PR_TIME("Matrix::copyTo");
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::copyTo | Error: Matrices must be the same size)");
#ifdef PR_BLIS
    bli_copym(&m_object, &(a_rhs.m_object));
#else
    // If we're copying a transpose, need to change order for mem linear copy
    MatrixStorageType not_storage = 
      (m_storage==ROW_MAJOR) ? COL_MAJOR : ROW_MAJOR;
    a_rhs.m_storage = (m_transpose) ? not_storage : m_storage;
    proto_memcpy<HOST, HOST>(m_data.get(), a_rhs.m_data.get(), linearSize());
#endif
}

template<typename T>
Matrix<T>
Matrix<T>::operator+(const Matrix<T>& a_rhs) const
{
    PR_TIME("matrix::operator+");
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
    Matrix<T> ret(M(), N());
#ifdef PR_BLIS
    this->copyTo(ret);
    bli_addm(&(a_rhs.m_object), &(ret.m_object));
#else
    for (int ii = 0; ii < size(); ii++)
    {
        ret.set(ii, this->get(ii) + a_rhs.get(ii));
    }
#endif
    return ret;
}

template<typename T>
void
Matrix<T>::operator+=(const Matrix<T>& a_rhs)
{
    PR_TIME("matrix::operator+=");
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
#ifdef PR_BLIS
    bli_addm(&(a_rhs.m_object), &m_object);
#else
    for (int ii = 0; ii < size(); ii++)
    {
        this->set(ii, this->get(ii) + a_rhs.get(ii));
    }
#endif
}

template<typename T>
Matrix<T>
Matrix<T>::operator+(T a_scale) const
{
    PR_TIME("matrix::operator+(scalar)");
    Matrix<T> ret(M(), N());
    ret.set(a_scale);
    ret += (*this);
    return ret;
}

template<typename T>
void
Matrix<T>::operator+=(T a_scale)
{
    PR_TIME("matrix::operator+=(scalar)");
    Matrix<T> tmp(M(), N());
    tmp.set(a_scale);
    (*this) += tmp;
}

template<typename T>
Matrix<T>
Matrix<T>::operator-(const Matrix<T>& a_rhs) const
{
    PR_TIME("matrix::operator-");
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
    Matrix<T> ret(M(), N());
#ifdef PR_BLIS
    this->copyTo(ret);
    bli_subm(&(a_rhs.m_object), &(ret.m_object));
#else
    for (int ii = 0; ii < size(); ii++)
    {
        ret.set(ii, this->get(ii) - a_rhs.get(ii));
    }

#endif
    return ret;
}

template<typename T>
void
Matrix<T>::operator-=(const Matrix<T>& a_rhs)
{
    PR_TIME("matrix::operator-=");
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
#ifdef PR_BLIS
    bli_subm(&(a_rhs.m_object), &m_object);
#else
    for (int ii = 0; ii < size(); ii++)
    {
        this->set(ii, this->get(ii) - a_rhs.get(ii));
    }
#endif
}

template<typename T>
Matrix<T>
Matrix<T>::operator-(T a_scale) const
{
    return (*this) + (-a_scale);
}

template<typename T>
void
Matrix<T>::operator-=(T a_scale)
{
    (*this) += (-a_scale);
}

template<typename T>
Matrix<T>
Matrix<T>::operator*(const Matrix<T>& a_rhs) const
{
    PR_TIME("Matrix::operator* (Matrix-Matrix multiply)");
    PROTO_ASSERT(N() == a_rhs.M(),
            "Matrix::operator* | Error: Incompatible Matrix dimensions.");
    Matrix<T> ret(M(), a_rhs.N());
#ifdef PR_BLIS
    if (a_rhs.isVector())
    {
        bli_gemv(&BLIS_ONE, &m_object, &(a_rhs.m_object), &BLIS_ZERO, &(ret.m_object));
    } else {
        bli_gemm(&BLIS_ONE, &m_object, &(a_rhs.m_object), &BLIS_ZERO, &(ret.m_object));
    }
#else
    
    PROTO_ASSERT(m_storage == COL_MAJOR,
            "Matrix multiply is not yet implemented for ROW_MAJOR storage");
    PROTO_ASSERT(a_rhs.m_storage == COL_MAJOR,
            "Matrix multiply is not yet implemented for ROW_MAJOR storage");

    int m1 =  M();
    int n1 =  N(); //== a_right.m_nrow
    int n2 =  a_rhs.N();
    int stride1 = m_numRows;
    int stride2 = a_rhs.m_numRows;

    char TRANS_L = 'N';
    if (m_transpose) { TRANS_L = 'T'; }
    char TRANS_R = 'N';
    if (a_rhs.m_transpose) { TRANS_R = 'T'; }
    double ALPHA = 1.;
    double BETA = 0.;
    PROTO_LAPACK(GEMM,gemm)(&TRANS_L, &TRANS_R, &m1, &n2, &n1, &ALPHA,
                           (T*)(*this).data(), &stride1,
                           (T*)a_rhs.data(), &stride2, &BETA,
                           ret.data(), &m1);
#endif
    return ret;
}

template<typename T>
Matrix<T>
Matrix<T>::operator*(T a_scale) const
{
    PR_TIME("Matrix::operator* (scalar)");
    Matrix<T> ret(M(), N());
    if (a_scale == 0)
    {
        ret.set(0);
    } else if (a_scale == 1)
    {
        this->copyTo(ret);
    } else {
#ifdef PR_BLIS
        obj_t alpha;
        bli_obj_create_1x1(blisDatatype(), &alpha);
        bli_setsc(a_scale, 0, &alpha);
        bli_scal2m(&alpha, &m_object, &(ret.m_object));
#else
        for (int ii = 0; ii < size(); ii++)
        {
            ret.set(ii, this->get(ii)*a_scale);
        }
#endif
    }
    return ret;
}

template<typename T>
void
Matrix<T>::operator*=(T a_scale)
{
    PR_TIME("Matrix::operator*= (scalar)");
    if (a_scale == 0)
    {
        this->set(0);
    } else if (a_scale == 1)
    {
        return;
    } else {
#ifdef PR_BLIS
        obj_t alpha;
        bli_obj_create_1x1(blisDatatype(), &alpha);
        bli_setsc(a_scale, 0, &alpha);
        bli_scalm(&alpha, &m_object);
#else
        for (int ii = 0; ii < size(); ii++)
        {
            this->set(ii, this->get(ii)*a_scale);
        }
#endif
    }
}

template<typename T>
T*
Matrix<T>::data()
{
#ifdef PR_BLIS
    return (T*)bli_obj_buffer(&m_object);
#else
    return m_data.get();
#endif
}

template<typename T>
const T*
Matrix<T>::data() const
{
#ifdef PR_BLIS
    return (T*)bli_obj_buffer(&m_object);
#else
    return m_data.get();
#endif
}

template<typename T>
const Matrix<T>
Matrix<T>::transpose() const
{
#ifdef PR_BLIS
    const Matrix<T> ret = slice(M(), N(), 0, 0);
    Matrix<T>* mut = const_cast<Matrix<T>*>(&ret);
    bli_obj_toggle_trans(&(mut->m_object));
    return ret;
#else
    Matrix<T> ret(M(), N());
    ret.m_data = m_data;
    ret.m_transpose = true;
    ret.m_alias = true;
    return ret;
#endif
}

template<typename T>
double
Matrix<T>::norm1() const
{
    PR_TIME("Matrix::norm1");
    double maxColSum = 0;
    for (int jj = 0; jj < N(); jj++)
    {
        double colSum = 0;
        for (int ii = 0; ii < M(); ii++)
        {
            colSum += std::abs(get(ii,jj));
        }
        maxColSum = std::max(colSum, maxColSum);
    }
    return maxColSum;
}
    
template<typename T>
double
Matrix<T>::normInf() const
{    
    PR_TIME("Matrix::normInf");
    double maxRowSum = 0;
    for (int ii = 0; ii < M(); ii++)
    {
        double rowSum = 0;
        for (int jj = 0; jj < N(); jj++)
        {
            rowSum += std::abs(get(ii,jj));
        }
        maxRowSum = std::max(rowSum, maxRowSum);
    }
    return maxRowSum;
}
template<typename T>
double
Matrix<T>::max() const
{
    PR_TIME("Matrix::max");
    double maxVal = get(0,0);
    for (int ii = 0; ii < M(); ii++)
    for (int jj = 0; jj < N(); jj++)
    {
        maxVal = std::max(maxVal, get(ii,jj));
    }
    return maxVal;
}

template<typename T>
double
Matrix<T>::absMax() const
{
    PR_TIME("Matrix::absMax");
    double maxVal = std::abs(get(0,0));
    for (int ii = 0; ii < M(); ii++)
    for (int jj = 0; jj < N(); jj++)
    {
        maxVal = std::max(maxVal, std::abs(get(ii,jj)));
    }
    return maxVal;
}

template<typename T>
double
Matrix<T>::min() const
{
    PR_TIME("Matrix::min");
    double minVal = get(0,0);
    for (int ii = 0; ii < M(); ii++)
    for (int jj = 0; jj < N(); jj++)
    {
        minVal = std::min(minVal, get(ii,jj));
    }
    return minVal;
}

template<typename T>
Matrix<T>
Matrix<T>::inverse() const
{
    if (M() == N())
    {
        PR_TIME("Matrix::inverse");
        Matrix<T> ret(M(), N());
        int n = M();
        int nn = n*n;
        const T* srcBuffer = data();
        T* dstBuffer = ret.data();
        proto_memcpy<HOST, HOST>(srcBuffer, dstBuffer, nn*sizeof(T));
        int* pivots = (int*)proto_malloc<HOST>(n*sizeof(int));
        T* work = (T*)proto_malloc<HOST>(nn*sizeof(T));
        int info;
        PROTO_LAPACK(GETRF,getrf)(&n, &n, dstBuffer, &n, pivots, &info);
        PROTO_ASSERT(info <= 0,
                "Matrix::inverse | Error: Matrix LU factorization failed. Matrix may be singular.");
        PROTO_LAPACK(GETRI,getri)(&n, dstBuffer, &n, pivots, work, &nn, &info);
        PROTO_ASSERT(info <= 0,
                "Matrix::inverse | Error: Matrix inversion failed. Matrix may be singular.");
        proto_free<HOST>(work);
        proto_free<HOST>(pivots);
        return ret;
    } else {
        return pseudoInverse(); 
    }
}

template<typename T>
Matrix<T>
Matrix<T>::pseudoInverse() const
{
    PR_TIME("Matrix::pseudoInverse");
    Matrix<T> U, S, VT;
    svd(U, S, VT);

    T* s = S.data();
    T* u = U.data();

    int m = U.M();
    int n = S.M();
    int incx = 1;
#pragma omp parallel for
    for(int i=0; i<n; i++)
    {
        double ss;
        if(s[i] > 1.0e-9) { ss=1.0/s[i]; }
        else { ss=s[i]; }
        PROTO_LAPACK(SCAL,scal)(&m, &ss, &u[i*m], &incx);
    }
   
    auto UT = U.transpose();
    auto V = VT.transpose();
    return V*UT;
}

template<typename T>
void
Matrix<T>::svd(Matrix<T>& a_U, Matrix<T>& a_S, Matrix<T>& a_VT) const
{
    PR_TIME("Matrix::SVD");
    char JOB = 'S';
    int m = M();
    int n = N();
    int LDA  = m;
    int LDU  = m;
    int LDVT = std::min(m,n);
    int INFO;

    a_U.define(LDU,LDVT);
    a_S.define(LDVT,1);
    a_VT.define(LDVT, n);

    // Create a copy of this (GESVD overwrites the input matrix)
    Matrix<T> A(M(), N());
    this->copyTo(A);

    // query and allocate optimal workspace 
    int LWORK = -1;
    double LWORK_OPT;
    PROTO_LAPACK(GESVD, gesvd)(&JOB, &JOB, &m, &n, A.data(), &LDA, a_S.data(),
            a_U.data(), &LDU, a_VT.data(), &LDVT, &LWORK_OPT, &LWORK, &INFO);
    
    LWORK = (int)LWORK_OPT;
    double* WORK = (T*)proto_malloc<HOST>(LWORK*sizeof(double));
    PROTO_LAPACK(GESVD, gesvd)(&JOB, &JOB, &m, &n, A.data(), &LDA, a_S.data(),
            a_U.data(), &LDU, a_VT.data(), &LDVT, WORK, &LWORK, &INFO);

    PROTO_ASSERT(INFO <= 0, "Matrix::svd | Error: SVD algorithm failed to converge.");
}

template<typename T>
unsigned int
Matrix<T>::size() const
{
    return M()*N();
}

template<typename T>
size_t
Matrix<T>::linearSize() const
{
    return size()*sizeof(T);
}

template<typename T>
bool 
Matrix<T>::isVector() const
{
    return ((numCols() == 1) || (numRows() == 1));
}

template<typename T>
bool 
Matrix<T>::isScalar() const
{
    return ((numCols() == 1) && (numRows() == 1));
}

template<typename T>
bool 
Matrix<T>::isTranspose() const
{
#ifdef PR_BLIS
    auto t = bli_obj_conjtrans_status(&m_object);
    if (t == BLIS_NO_TRANSPOSE)
    {
        return false;
    } else if (t == BLIS_TRANSPOSE)
    {
        return true;
    }
#else
    return m_transpose;
#endif
}

template<typename T>
void 
Matrix<T>::printBuffer(std::string a_format) const
{
    const T* buffer = this->data();
    for (int ii = 0; ii < size(); ii++)
    {
        Proto::pout() << buffer[ii] << ", ";
    }
    Proto::pout() << std::endl;
}

template<typename T>
void 
Matrix<T>::print(std::string a_format) const
{
#ifdef PR_BLIS
    Matrix<T> A(M(), N());
    this->copyTo(A); //this make sure transposed matrices print correctly
    if (a_format == "")
    {
        switch (blisDatatype())
        {
            case BLIS_DOUBLE:
                bli_printm("", &(A.m_object), "%4.2f", "");
                break;
            default:
                bli_printm("", &(A.m_object), "%i", "");
                break;

        }
    } else {
        bli_printm("", &m_object, a_format.c_str(), "");
    }
#else
    for (int ii = 0; ii < M(); ii++)
    {
        Proto::pout() << std::endl;
        for (int jj = 0; jj < N(); jj++)
        {
            Proto::pout() << get(ii,jj);
            if (jj != N()-1) {Proto::pout() << ", "; }
        }
    }
    Proto::pout() << std::endl;
#endif
}

template<typename T>
double
Matrix<T>::conditionNumber() const
{
    auto Inv = inverse();
    return this->normInf()*Inv.normInf();
}

template<typename T>
void
Matrix<T>::checkConditionNumber() const
{
    double inverse = 1.0/conditionNumber();
    double small = 1.0e-6;
    double reallysmall = 1.0e-15;
    if(inverse < reallysmall)
    {
        Proto::pout() << "matrix is poorly conditioned: 1/condition_number = " << inverse << endl;
    }
    else if(inverse < small)
    {
        Proto::pout() << "matrix is marginally conditioned: 1/condition_number = " << inverse << endl;
    }
    else
    {
        Proto::pout() << "matrix might be OK: 1/condition_number = " << inverse << endl;
    }
}

template<typename T>
void solveLS(Matrix<T>& A, Matrix<T>& B)
{
    PR_TIME("solveLeastSquares");
    // TODO - check that the sizes of A, B and C are compatible
    int M = A.M();
    int N = A.N();
    int NRHS = B.N();
    int LDA = M;
    int LDB = std::max(M,N);
    PROTO_ASSERT(B.M() == M,
            "solveLS | Error: Incompatible matrix dimensions. A.M() == %u, b.M() = %u.",
            A.M(), B.M());

    int LWORK[2] = {1,1};
    LWORK[0] = 2*M*N;

    Matrix<T> WORK(2*M*N, 1);
    WORK.set(0.0);

    char TRANS = 'N';
    if (A.isTranspose()) { TRANS = 'T'; }
    int INFO;
   
    PROTO_LAPACK(GELS,gels)(&TRANS, &M, &N, &NRHS, A.data(), &LDA, 
            B.data(), &LDB, WORK.data(), LWORK, &INFO);

    PROTO_ASSERT(INFO == 0,
            "solveLS | Error: GELS returned an error flag. Matrix may be singular.");
}

template<typename T>
void solveRRLS(Matrix<T>& A, Matrix<T>& b)
{
    // - check that the sizes of A, b are compatible
    int M = A.M();
    int N = A.N();
    int NRHS =  b.N();
    PROTO_ASSERT(b.M() == M,
            "solveRRLS | Error: Incompatible matrix domensions. A.M() = %u, b.M() = %u.",
            A.M(), b.M());
    PROTO_ASSERT(M >= N, "solveRRLS | Error: System is not over-determined.");

    // - Compute A transpose
    auto At = A.transpose(); 

    // - Calculate QR factorization of A' with DGEQRF
    int LWORK[2] = {1,1};
    LWORK[0] = 2*M*N;
    Matrix<T> WORK(LWORK[0], 1);
    WORK.setVal(0.0);

    Matrix<T> TAU(N, 1);
    TAU.setVal(0.0);

    int INFO;
    PROTO_LAPACK(GEQRF,geqrf)(&N, &M, At.dataPtr(), &N,
            TAU.dataPtr(), WORK.dataPtr(), LWORK, &INFO);
    
    PROTO_ASSERT(INFO == 0,
            "solveRRLS | Error: GEQRF returned an error flag. Matrix may be singular.");

    char SIDE = 'L';
    char TRANS = 'T';
    PROTO_LAPACK(ORMQR,ormqr)(&SIDE, &TRANS, &N, &NRHS, &N, 
            At.data(), &N, TAU.data(), b.data(), &N,
            WORK.data(), LWORK, &INFO);

    PROTO_ASSERT(INFO == 0,
            "solveRRLS | Error: ORMQR returned an error flag. Matrix may be singular.");

    // - Solve R x = (Q' * b) using DTRTRS
    char UPLO = 'U';
    TRANS = 'N';
    char DIAG = 'N';
    PROTO_LAPACK(TRTRS,trtrs)(&UPLO, &TRANS, &DIAG, &N, &NRHS, 
            At.data(), &N, b.data(), &N, &INFO);
    
    PROTO_ASSERT(INFO == 0,
            "solveRRLS | Error: TRTRS returned an error flag. Matrix may be singular.");
}
