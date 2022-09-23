
template<typename T>
Matrix<T>::Matrix(unsigned int a_numRows, unsigned int a_numCols)
{
    bli_obj_create(blisDatatype(), a_numRows, a_numCols, 0, 0, &m_object);
    m_alias = false;
}

template<typename T>
Matrix<T>::Matrix(std::shared_ptr<T>& a_buffer, unsigned int a_numRows, unsigned int a_numCols)
{
    bli_obj_create_without_buffer(blisDatatype(), a_numRows, a_numCols, &m_object);
    bli_obj_attach_buffer(a_buffer.get(), 1, a_numRows, 0, &m_object);
    m_alias = true;
}

template<typename T>
Matrix<T>::Matrix(const Matrix<T>&& a_matrix)
{
    m_alias = a_matrix.m_alias;
    m_object = std::move(a_matrix.m_object);
}

template<typename T>
Matrix<T>::~Matrix()
{
    if (!m_alias)
    {
        bli_obj_free(&m_object);
    }
}

template<typename T>
Matrix<T>&
Matrix<T>::operator=(const Matrix<T>&& a_matrix)
{
    m_alias = a_matrix.m_alias;
    m_object = std::move(a_matrix.m_object);
    return *this;
}

template<typename T>
unsigned int
Matrix<T>::numRows() const
{
    return bli_obj_length(&m_object);
}

template<typename T>
unsigned int
Matrix<T>::numCols() const
{
    return bli_obj_width(&m_object);
}

template<typename T>
Matrix<T>
Matrix<T>::slice(
        unsigned int a_numRows, unsigned int a_numCols,
        unsigned int a_rowOffset, unsigned int a_colOffset)
{
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
Matrix<T>::row(unsigned int a_rowIndex) const
{
}

template<typename T>
Matrix<T>
Matrix<T>::col(unsigned int a_colIndex) const
{
}

