
template<typename T>
Matrix<T>::Matrix(unsigned int a_numRows, unsigned int a_numCols)
{
    PR_TIME("Matrix(numRows, numCols)");
    bli_obj_create(blisDatatype(), a_numRows, a_numCols, 0, 0, &m_object);
    m_alias = false;
}

template<typename T>
Matrix<T>::Matrix(std::shared_ptr<T>& a_buffer, unsigned int a_numRows, unsigned int a_numCols)
{
    PR_TIME("Matrix(buffer, numRows, numCols)");
    bli_obj_create_without_buffer(blisDatatype(), a_numRows, a_numCols, &m_object);
    bli_obj_attach_buffer(a_buffer.get(), 1, a_numRows, 0, &m_object);
    m_alias = true;
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>&& a_matrix)
{
    PR_TIME("Matrix(Matrix&&)");
    m_alias = a_matrix.m_alias;
    m_object = std::move(a_matrix.m_object);
}

template<typename T>
Matrix<T>::~Matrix()
{
    PR_TIME("~Matrix()");
    if (!m_alias)
    {
        bli_obj_free(&m_object);
    }
}

template<typename T>
Matrix<T>&
Matrix<T>::operator=(const Matrix<T>&& a_matrix)
{
    PR_TIME("Matrix::operator=(Matrix&&)");
    m_alias = a_matrix.m_alias;
    m_object = std::move(a_matrix.m_object);
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
    T real;
    T imag; //unused
    //TODO: make this work properly if T is complex
    bli_getijm(a_row, a_col, &m_object, &real, &imag);
    return real;
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
    bli_setijm(a_val, 0, a_row, a_col, &m_object);
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
}

template<typename T>
unsigned int
Matrix<T>::numRows() const
{
    PR_TIME("Matrix::numRows");
    return bli_obj_length(&m_object);
}

template<typename T>
unsigned int
Matrix<T>::numCols() const
{
    PR_TIME("Matrix::numCols");
    return bli_obj_width(&m_object);
}

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

template<typename T>
void
Matrix<T>::copyTo(Matrix<T>& a_rhs) const
{
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::copyTo | Error: Matrices must be the same size)");
    bli_copym(&m_object, &(a_rhs.m_object));
}

template<typename T>
Matrix<T>
Matrix<T>::operator+(const Matrix<T>& a_rhs) const
{
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
    Matrix<T> ret(M(), N());
    this->copyTo(ret);
    bli_addm(&(ret.m_object), &(a_rhs.m_object));
    return ret;
}

template<typename T>
void
Matrix<T>::operator+=(const Matrix<T>& a_rhs)
{
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
    bli_addm(&m_object, &(a_rhs.m_object));
}

template<typename T>
Matrix<T>
Matrix<T>::operator+(T a_scale) const
{
    Matrix<T> ret(M(), N());
    ret.set(a_scale);
    ret += (*this);
    return ret;
}

template<typename T>
void
Matrix<T>::operator+=(T a_scale)
{
    Matrix<T> tmp(M(), N());
    tmp.set(a_scale);
    (*this) += tmp;
}

template<typename T>
Matrix<T>
Matrix<T>::operator-(const Matrix<T>& a_rhs) const
{
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
    Matrix<T> ret(M(), N());
    this->copyTo(ret);
    bli_subm(&(ret.m_object), &(a_rhs.m_object));
    return ret;
}

template<typename T>
void
Matrix<T>::operator-=(const Matrix<T>& a_rhs)
{
    PROTO_ASSERT(this->M() == a_rhs.M() && this->N() == a_rhs.N(),
            "Matrix::operator+ | Error: Matrices must be the same size)");
    bli_subm(&m_object, &(a_rhs.m_object));
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
    PROTO_ASSERT(N() == a_rhs.M(),
            "Matrix::operator* | Error: Incompatible Matrix dimensions.");
    Matrix<T> ret(M(), a_rhs.N());
    if (a_rhs.isVector())
    {
        bli_gemv(&BLIS_ONE, &m_object, &(a_rhs.m_object), &BLIS_ZERO, &(ret.m_object));
    } else {
        bli_gemm(&BLIS_ONE, &m_object, &(a_rhs.m_object), &BLIS_ZERO, &(ret.m_object));
    }
    return ret;
}

template<typename T>
void
Matrix<T>::operator*=(const Matrix<T>& a_rhs)
{
    (*this) = (*this)*a_rhs;
}

template<typename T>
Matrix<T>
Matrix<T>::operator*(T a_scale) const
{
    Matrix<T> ret(M(), N());
    if (a_scale == 0)
    {
        ret.set(0);
    } else if (a_scale == 1)
    {
        this->copyTo(ret);
    } else {
        obj_t alpha;
        bli_obj_create_1x1(blisDatatype(), &alpha);
        bli_setsc(a_scale, 0, &alpha);
        bli_scal2m(&alpha, &m_object, &(ret.m_object));
    }
    return ret;
}

template<typename T>
void
Matrix<T>::operator*=(T a_scale)
{
    if (a_scale == 0)
    {
        this->set(0);
    } else if (a_scale == 1)
    {
        return;
    } else {
        obj_t alpha;
        bli_obj_create_1x1(blisDatatype(), &alpha);
        bli_setsc(a_scale, 0, &alpha);
        bli_scalm(&alpha, &m_object);
    }
    return ret;
}

template<typename T>
Matrix<T>
Matrix<T>::transpose() const
{
    Matrix<T> ret = slice(M(), N(), 0, 0);
    bli_obj_toggle_trans(&(ret.m_object));
    return ret;
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
void 
Matrix<T>::print(std::string a_format) const
{
    if (a_format == "")
    {
        switch (blisDatatype())
        {
            case BLIS_DOUBLE:
                bli_printm("", &m_object, "%4.1f", "");
                break;
            default:
                bli_printm("", &m_object, "%i", "");
                break;

        }
    } else {
        bli_printm("", &m_object, a_format.c_str(), "");
    }
}
