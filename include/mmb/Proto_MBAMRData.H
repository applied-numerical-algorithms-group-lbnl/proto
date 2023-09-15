#pragma once
#ifndef _PROTO_MB_AMR_DATA_
#define _PROTO_MB_AMR_DATA_

#include "Proto_MBAMRGrid.H"
#include "Proto_MBLevelBoxData.H"

namespace Proto
{
    /// Multiblock AMR Data
    template<typename T, unsigned int C, MemType MEM, Centering CTR=PR_CELL>
    class MBAMRData
    {
        public:
        
        inline MBAMRData();
        inline MBAMRData(
                const MBAMRGrid& a_grid,
                Array<Point, DIM+1> a_ghost,
                Point a_boundGhost = Point::Zeros());
        
        inline void define(
                const MBAMRGrid& a_grid,
                Array<Point, DIM+1> a_ghost,
                Point a_boundGhost = Point::Zeros());
        
        inline void setVal(T a_value);
        inline void setRandom(T a_low, T a_high);
        inline void exchange();
        inline Array<Point, DIM+1> ghost() const { return m_ghost; }
        inline MBLevelBoxData<T,C,MEM,CTR>& getLevel(unsigned int a_level);
        inline const MBLevelBoxData<T,C,MEM,CTR>& getLevel(unsigned int a_level) const;
        inline AMRData<T,C,MEM,CTR>& getBlock(unsigned int a_block);
        inline const AMRData<T,C,MEM,CTR>& getBlock(unsigned int a_block) const;
        inline const MBAMRGrid& grid() const { return m_grid; } 
        inline int numLevels() const {return m_grid.numLevels(); }
        inline int numBlocks() const {return m_grid.numBlocks(); }

        // Operators
        inline void averageDown();        
        inline void increment(const MBAMRData<T,C,MEM,CTR>& a_data, T a_scale);
        inline void operator*=(T a_value);
        inline void operator+=(T a_value);
        inline T integrate(Array<T,DIM> a_cdx, unsigned int a_c = 0);
        inline T integrate(T a_cdx, unsigned int a_c = 0);
        inline T integrateAbs(Array<T,DIM> a_cdx, unsigned int a_c = 0);
        inline T integrateAbs(T a_cdx, unsigned int a_c = 0);
        inline T absMax(unsigned int a_c = 0);
        inline void copyTo(MBAMRData& a_rhs);
        
        private:
        
        Array<Point, DIM+1> m_ghost;
        const MBAMRGrid*    m_grid;
        std::vector<std::shared_ptr<MBLevelBoxData<T, C, MEM, CTR>>>    m_levelData;
        std::vector<std::shared_ptr<AMRData<T, C, MEM, CTR>>>           m_amrData;
    };

#include "implem/Proto_MBAMRDataImplem.H"
} // end namespace Proto
#endif //end include guard
