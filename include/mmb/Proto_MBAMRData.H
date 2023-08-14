#pragma once
#ifndef _PROTO_MB_AMR_DATA_
#define _PROTO_MB_AMR_DATA_

#include "Proto_MBLevelBoxData.H"

namespace Proto
{
    /// Multiblock AMR Data
    template<typename T, unsigned int C, MemType MEM, Centering CTR=PR_CELL>
    class MBAMRData
    {
        public:


        private:
        
        Array<Point, DIM+1>     m_ghost;
        const MBAMRLayout*      m_layout;
        std::vector<std::shared_ptr<MBLevelBoxData<T, C, MEM, CTR>>> m_levelData;
        std::vector<std::shared_ptr<AMRData<T, C, MEM, CTR>>> m_amrData;
    };

#include "implem/Proto_MBAMRDataImplem.H"
} // end namespace Proto
#endif //end include guard
