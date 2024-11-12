#pragma once
#ifndef _PROTO_REGISTER_H_
#define _PROTO_REGISTER_H_
#include "Proto.H"
#include "Proto_MemType.H"
namespace Proto
{
    template<typename T, unsigned int C, MemType MEMTYPE>
    struct Register
    {
        std::shared_ptr<BoxData<T,C,MEMTYPE> > m_data;
        unsigned int m_dir;
        Side::LoHiSide m_side;
        //unsigned int m_lohi;
        Register(
                std::shared_ptr<BoxData<T,C,MEMTYPE> > a_data,
                unsigned int a_dir,
                Side::LoHiSide a_side)
        {
            m_data = a_data;
            m_dir = a_dir;
            m_side = a_side;
        };

        Register(
                std::shared_ptr<BoxData<T,C,MEMTYPE>> a_data,
                Point a_dir)
        {
            PROTO_ASSERT(a_dir.codim() == 1,
                "Register constructor | Error: Direction Point must have codimension 1");
            m_data = a_data;
            for (int dir = 0; dir < DIM; dir++)
            {
                if (a_dir[dir] > 0)
                {
                    m_dir = dir;
                    m_side = Side::Hi;
                } else if (a_dir[dir] < 0)
                {
                    m_dir = dir;
                    m_side = Side::Lo;
                } else {
                    PROTO_ASSERT(false,
                            "Register constructor | Error: Data corruption.");
                }
            }
        }
    };
}
#endif
