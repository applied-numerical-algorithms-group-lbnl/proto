
#pragma once
#ifndef _PROTO_CELL_AVERAGE_MATRIX_PRODUCT_OP_
#define _PROTO_CELL_AVERAGE_MATRIX_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    /// @brief Computes the cell averaged product of two matrix-valued BoxData to 4th order
    /// @tparam T 
    template <class T>
    class CellAverageMatrixProductOp : public MatrixProductOp<T, CellAverageProductOp>
    {

    private:

    public:
        inline CellAverageMatrixProductOp()
        {
            this->_productOp = std::make_shared<CellAverageProductOp<T>>();
        }
    };
}
#endif