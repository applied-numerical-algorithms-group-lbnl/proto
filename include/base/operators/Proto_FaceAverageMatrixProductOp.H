
#pragma once
#ifndef _PROTO_FACE_AVERAGE_MATRIX_PRODUCT_OP_
#define _PROTO_FACE_AVERAGE_MATRIX_PRODUCT_OP_

#include "Proto.H"

namespace Proto
{
    template <class T>
    class FaceAverageMatrixProductOp : public MatrixProductOp<T, FaceAverageProductOp>
    {

    private:
        int _normDir;

    public:
        inline FaceAverageMatrixProductOp(int normDir)
        {
            _normDir = normDir;
            this->_productOp = std::make_shared<FaceAverageProductOp<T>>(normDir);
        }
    };
}
#endif