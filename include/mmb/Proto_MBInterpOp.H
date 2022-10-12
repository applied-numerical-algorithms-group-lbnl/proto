#pragma once
#ifndef _PROTO_MB_INTERP_OP_
#define _PROTO_MB_INTERP_OP_

#include "Proto.H"
#include "Proto_MBLevelBoxData.H" //for MBDataPoint definition

namespace Proto
{
    template<typename T>
    class MBInterpOp
    {
        public:

        template<unsigned int C, MemType MEM>
        inline std::array<T, C> execute(
                MBLevelBoxData<T, C, MEM>& a_dst,
                MBLevelBoxData<T, C, MEM>& a_src);
        
        inline unsigned int size() const;

        private:

        MBDataPoint m_dst;
        std::vector<MBDataPoint> m_srcs;
        std::vector<T> m_coefs;
    }
#include "implem/Proto_MBInterpOpImplem.H"
} //end namespace Proto
#endif // end include guard
