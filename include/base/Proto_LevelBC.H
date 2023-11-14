#pragma once
#ifndef __PROTO_LEVEL_BC__
#define __PROTO_LEVEL_BC__

#include "Proto_LevelBoxData.H"
namespace Proto {
    
#ifdef PR_MMB
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
        class MBLevelBoxData;
#endif

    /// Storage for the LevelBCOps on each face 
    template<typename T, unsigned int C, MemType MEM=MEMTYPE_DEFAULT, Centering CTR=PR_CELL>
    class LevelBC {
        public:

        // Constructor
        inline LevelBC(){};
        inline LevelBC(DisjointBoxLayout& a_layout);
        inline void define(DisjointBoxLayout& a_layout);
        
        // Basic interface
        inline virtual void init() {}
        inline virtual void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const {}
        #ifdef PR_MMB
        inline virtual void apply(MBLevelBoxData<T,C,MEM,CTR>& a_data) const {}
        #endif

        // Helper methods (not yet implemented)
        inline const DisjointBoxLayout& layout() { return m_layout; }
        inline void setFlux(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_fluxValue, int a_comp, int a_coordinate, Side::LoHiSide a_side);
        inline void setValue(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_fluxValue, int a_comp, int a_coordinate, Side::LoHiSide a_side);
        
        // ... many useful things here ...

        inline void setTime(T a_time) const;
        inline T time() const { return m_time; }

        // Scalar constant interface
        /// Set Number of Varying Constants
        inline void setNumConstants(int a_numConst) const;
        /// Set Specific Varying Constant
        inline void setConstant(int a_index, T a_value) const;
        /// Get Specific Varying Constant
        inline T getConstant(int a_index) const;

        // Fixed data interface
        /// Set Number of Variable LevelData
        inline void setNumLevelData(int a_numFixed) const;
        /// Set Variable LevelData
        inline void setLevelData(int a_index, const LevelBoxData<T,C,MEM,CTR>& a_data) const;
        /// Get Variable LevelData
        inline const LevelBoxData<T,C,MEM,CTR>& getLevelData(int a_index) const;

        // RK interface
        inline void setRKStage(int a_stage) const { m_rkStage = a_stage; }
        inline void setNumRKStages(int a_numStages) const;
        inline void setRKStageData(int a_stage, const LevelBoxData<T,C,MEM,CTR>& a_data) const;
        inline const LevelBoxData<T,C,MEM,CTR>& getRKStageData(int a_index) const;
        inline int rkStage() { return m_rkStage; }

        private:

        DisjointBoxLayout m_layout;

        mutable T m_time;
        mutable int m_rkStage;
        mutable std::vector<T> m_constData;
        mutable std::vector<const LevelBoxData<T,C,MEM,CTR>*> m_fixedData;
        mutable std::vector<const LevelBoxData<T,C,MEM,CTR>*> m_rkStageData;
    };

    #include "implem/Proto_LevelBCImplem.H"
} // end namespace Proto
#include "Proto_LevelBCLib.H" // intentionally outside the namespace
#endif //end include guard
