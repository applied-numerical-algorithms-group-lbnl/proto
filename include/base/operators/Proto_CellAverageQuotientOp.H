
#pragma once
#ifndef _PROTO_CELL_AVERAGE_QUOTIENT_OP_
#define _PROTO_CELL_AVERAGE_QUOTIENT_OP_

#include "Proto.H"

namespace Proto
{
    /// @brief Computes the cell-averaged quotient of two cell-averaged BoxData to 4th order
    /// @tparam T 
    template <class T>
    class CellAverageQuotientOp : public QuotientAverageOp<T>
    {

    private:
        Array<Stencil<T>, DIM> _diffs;

    public:
        inline CellAverageQuotientOp()
        {
            for (int dd = 0; dd < DIM; dd++)
            {
                Stencil<T> diff = 0.5*Shift::Basis(dd, 1) - 0.5*Shift::Basis(dd, -1);
                this->_diffOps.push_back(diff);
            }
        }
    };
}
#endif