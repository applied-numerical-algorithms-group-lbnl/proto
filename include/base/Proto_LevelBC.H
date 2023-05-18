#pragma once
#ifndef __PROTO_LEVEL_BC__
#define __PROTO_LEVEL_BC__

#include "Proto_LevelBoxData.H"

namespace Proto {
    
    /// Storage for the LevelBCOps on each face 
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class LevelBC {
        public:
        
        // Basic interface
        inline virtual void apply(LevelBoxData<T,C,MEM,CTR>& a_data) const {}

        // Helper methods (not yet implemented)
        inline void setFlux(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_fluxValue, int a_comp, int a_coordinate, Side::LoHiSide a_side);
        inline void setValue(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_fluxValue, int a_comp, int a_coordinate, Side::LoHiSide a_side);
        
        // ... many useful things here ...

        // Fixed data interface
        inline void setNumFixedData(int a_numFixed) const;
        inline void setFixedData(int a_index, const LevelBoxData<T,C,MEM,CTR>& a_data) const;
        inline const LevelBoxData<T,C,MEM,CTR>& getFixedData(int a_index) const;

        // RK interface
        inline void setRKStage(int a_stage) const { m_rkStage = a_stage; }
        inline void setNumRKStages(int a_numStages) const;
        inline void setRKStageData(int a_stage, const LevelBoxData<T,C,MEM,CTR>& a_data) const;
        inline const LevelBoxData<T,C,MEM,CTR>& getRKStageData(int a_index) const;
        inline int rkStage() { return m_rkStage; }

        private:
        mutable int m_rkStage;
        mutable std::vector<const LevelBoxData<T,C,MEM,CTR>*> m_fixedData;
        mutable std::vector<const LevelBoxData<T,C,MEM,CTR>*> m_rkStageData;
    };

    #include "implem/Proto_LevelBCImplem.H"
} // end namespace Proto
#include "Proto_LevelBCLib.H" // intentionally outside the namespace
#endif //end include guard
