#pragma once
#ifndef __PROTO_LEVEL_BC__
#define __PROTO_LEVEL_BC__

#include "Proto_LevelBoxData.H"
namespace Proto {
    
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

        
        // NOT IMPLEMENTED YET
        inline void fillGhost(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_value, int a_comp, int a_coordinate, Side::LoHiSide a_size);
        // NOT IMPLEMENTED YET
        inline void setFlux(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_fluxValue, int a_comp, int a_coordinate, Side::LoHiSide a_side);
        // NOT IMPLEMENTED
        inline void setValue(LevelBoxData<T,C,MEM,CTR>& a_data,
                T a_value, int a_comp, int a_coordinate, Side::LoHiSide a_side);
        
        inline T time() const { return m_time; }
        inline int rkStage() { return m_rkStage; }
        inline const DisjointBoxLayout& layout() { return m_layout; }

        // RK interface
        inline void setTime(T a_time) const;
        inline void setRKStage(int a_stage) const { m_rkStage = a_stage; }
        inline void setNumRKStages(int a_numStages) const;
        inline void setRKStageData(int a_stage, const LevelBoxData<T,C,MEM,CTR>& a_data) const;
        inline const LevelBoxData<T,C,MEM,CTR>& getRKStageData(int a_index) const;

        private:

        DisjointBoxLayout m_layout;
        mutable T m_time;
        mutable int m_rkStage;
        mutable std::vector<const LevelBoxData<T,C,MEM,CTR>*> m_rkStageData;
    };

    #include "implem/Proto_LevelBCImplem.H"
} // end namespace Proto
#include "Proto_LevelBCLib.H" // intentionally outside the namespace
#endif //end include guard
