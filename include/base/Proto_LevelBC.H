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
        inline virtual void apply(LevelBoxData<T,C,MEM>& a_data);

        // Fixed data interface
        inline void setNumFixedData(int a_numFixed);
        inline void setFixedData(int a_index, const LevelBoxData<T,C,MEM>& a_data);
        inline const LevelBoxData<T,C,MEM>& getFixedData(int a_index);

        // RK interface
        inline void setNumRKStages(int a_numStages);
        inline void setRKStageData(int a_stage, const LevelBoxData<T,C,MEM>& a_data);
        inline const LevelBoxData<T,C,MEM>& getRKStageData(int a_index);
        
        private:
        std::vector<std::weak_ptr<const LevelBoxData<T,C,MEM>>> m_fixedData;
        std::vector<std::weak_ptr<const LevelBoxData<T,C,MEM>>> m_rkStageData;

    };

    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class PeriodicBC : LevelBC
    {
        inline void apply(LevelBoxData<T,C,MEM>& a_data)
        {
            a_data.exchange();
        }
    }

    #include "implem/Proto_LevelBCImplem.H"
} // end namespace Proto
#endif //end include guard
