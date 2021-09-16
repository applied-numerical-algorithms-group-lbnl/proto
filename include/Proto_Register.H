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
        shared_ptr<BoxData<T,C,MEMTYPE> > m_data;
        unsigned int m_dir;
        unsigned int m_side;
        Register(
                shared_ptr<BoxData<T,C,MEMTYPE> > a_data,
                unsigned int a_dir,
                Side a_side)
        {
            m_data = a_data;
            m_dir = a_dir;
            m_side = a_side;
        };
    };
}
#endif
