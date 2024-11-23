#pragma once
#ifndef _PROTO_REGISTER_H_
#define _PROTO_REGISTER_H_
#include "Proto.H"
#include "Proto_MemType.H"
namespace Proto
{
    template<typename T, unsigned int C, MemType MEM>
    struct Register
    {
        std::shared_ptr<BoxData<T,C,MEM> > m_data;
        unsigned int m_dir;
        Side::LoHiSide m_side;
        Register(
                std::shared_ptr<BoxData<T,C,MEM> > a_data,
                unsigned int a_dir,
                Side::LoHiSide a_side)
        {
            m_data = a_data;
            m_dir = a_dir;
            m_side = a_side;
        };

        Register(
                std::shared_ptr<BoxData<T,C,MEM>> a_data,
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
                }
            }
        }

        BoxData<T,C,MEM>& data() { return *m_data; } 

        Point dir() const { return Point::Basis(m_dir, m_side); }

    };
}
#endif
