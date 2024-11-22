#pragma once
#ifndef _PROTO_MB_FLUX_REGISTER_
#define _PROTO_MB_FLUX_REGISTER_

#include "Proto.H"

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM>
    class MBLevelFluxRegisterTester;
    template<typename T, unsigned int C, MemType MEM>
    class MBLevelFluxRegister;
    template<typename T, unsigned int C, MemType MEM>
    class MBFluxRegisterCopier;

    template <typename T, unsigned int C, MemType MEM=MEMTYPE_DEFAULT>
    class MBLevelFluxRegisterTester
    {
    public:
        MBLevelFluxRegisterTester(const MBLevelFluxRegister<T, C, MEM> &a_register) : m_register(a_register) {}
        inline const std::vector<Register<T, C, MEM>>& getCoarseRegistersAtIndex(MBIndex a_coarseIndex) const;
        inline const std::vector<Register<T, C, MEM>>& getFineRegistersAtIndex(MBIndex a_fineIndex) const;

    private:
        const MBLevelFluxRegister<T, C, MEM> &m_register;
    };

    template<typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT> 
    class MBFluxRegisterCopierOp
    {
        friend class MBFluxRegisterCopier<T,C,MEM>;
        private:
        MBLevelFluxRegister<T, C, MEM>* m_register;

        public:
        inline MBFluxRegisterCopierOp(){};
        inline MBFluxRegisterCopierOp(MBLevelFluxRegister<T, C, MEM>& a_register);

        /// Get Serialized Size
        inline int linearSize(const Box& a_bx,
                const MBIndex& a_index) const;

        /// Serial Write To Buffer
        inline void linearOut(void* a_buf, const MBMotionItem& a_info) const;

        /// Serial Read From Buffer
        inline void linearIn(void* a_buf, const MBMotionItem& a_info) const;

        /// Local Copy
        inline void localCopy(const MBMotionItem& a_info) const;

    }; // end class MBFluxRegisterCopierOp

    template<typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT>
    class MBFluxRegisterCopier
    : public Copier<MBFluxRegisterCopierOp<T, C, MEM>, MBBoxPartition, MBBoxPartition, MEM, MEM>
    {
        public:

        inline MBFluxRegisterCopier(){};

        /// Build Copier Motion Plan
        inline void buildMotionPlans(MBFluxRegisterCopierOp<T, C, MEM>& a_op);
    }; //end class MBFluxRegisterCopier

    template<typename T, unsigned int C, MemType MEM=MEMTYPE_DEFAULT>
    class MBLevelFluxRegister
    {
        friend class MBLevelFluxRegisterTester<T,C,MEM>;
        friend class MBFluxRegisterCopierOp<T,C,MEM>;
        friend class MBFluxRegisterCopier<T,C,MEM>;

        public:
        MBLevelFluxRegister(){};
        MBLevelFluxRegister(const MBDisjointBoxLayout& a_dblCoarse,
                const MBDisjointBoxLayout&  a_dblFine,
                Point&          a_refRatio,
                Array<T, DIM>&  a_dxCoarse);
        
        MBLevelFluxRegister(const MBDisjointBoxLayout& a_dblCoarse,
                const MBDisjointBoxLayout&  a_dblFine,
                const std::vector<Point>&   a_refRatios,
                std::vector<Array<T, DIM>>& a_dxCoarse);

        void define(const MBDisjointBoxLayout& a_dblCoarse,
                const MBDisjointBoxLayout&  a_dblFine,
                Point&          a_refRatios,
                Array<T, DIM>&  a_dxCoarse);
        
        void define(const MBDisjointBoxLayout& a_dblCoarse,
                const MBDisjointBoxLayout&  a_dblFine,
                const std::vector<Point>&   a_refRatios,
                std::vector<Array<T, DIM>>& a_dxCoarse);

        void clearInternalData();
        void initCoarseRegisters(T a_val = 0);
        void initFineRegisters(T a_val = 0);
        void initIncrement(T a_val = 0);

        /// Increment Coarse
        /**
          Increments all the registers contained in the coarse patch a_di, direction a_dir. 
          Assumes that a_flux is defined on all the a_dir faces in the patch.
         */
        void incrementCoarseRegister(
                const BoxData<T,C,MEM>& a_flux,
                const MBIndex& a_di,
                const T& a_weight,
                unsigned int a_dir);

        /// Increment Fine
        /**
          Increments all the registers contained in the fine patch a_di, direction a_dir. 
          Assumes that a_flux is defined on all the a_dir faces in the patch.
         */
        void incrementFineRegister(const BoxData<T,C,MEM>& a_flux,
                const MBIndex& a_di,
                const T& a_weight,
                unsigned int a_dir); 

        /// Reflux
        /**
          Performs the reflux correction on input coarse data. 
         */
        void applyRefluxCorrection(MBLevelBoxData<T,C,MEM>& a_coarseData, const T& a_weight);
        

        
        /// compute the hash tileIndex for a coarse register.
        int tileIndex(const Box& a_bx,const MBIndex& a_di);

        void print() const;
      
        /// Access to layouts.
        MBDisjointBoxLayout crseLayout() const {return m_crseLayout;};
        MBDisjointBoxLayout fineLayout() const {return m_fineLayout;};
        
        protected:
        void validateInputs(const MBDisjointBoxLayout &a_dblCoarse,
                const MBDisjointBoxLayout &a_dblFine,
                const std::vector<Point> &a_refRatios,
                std::vector<Array<T, DIM>> &a_dxCoarse) const;
        void defineCoarse();
        void defineFine();
        void defineAverageStencils();
        

        BoxData<T, C, MEM> &sourceData(const Box &a_bx, const MBIndex &di);
        BoxData<T, C, MEM> &destData(const Box &a_bx, const MBIndex &di);

        MBDisjointBoxLayout m_crseLayout;
        MBDisjointBoxLayout m_fineLayout;
        std::vector<Array<T, DIM>> m_dxCrse;
        MBFluxRegisterCopier<T, C, MEM> m_copier;
        std::vector<Point> m_refRatios;

        // Registers for coarse data
        std::vector<std::vector<Register<T, C, MEM>>> m_crseRegisters;
        // Registers for fine data on coarse layout
        std::vector<std::vector<shared_ptr<BoxData<T, C, MEM>>>> m_crseIncrement;
        // Registers for fine data
        std::vector<std::vector<Register<T, C, MEM>>> m_fineRegisters;
        // crseIndices is a map because there are in general multiple registers for each (dir, side)
        std::vector<unordered_map<unsigned int, unsigned int>> m_crseIndices;
        // fineIndices are not a map because there is at most one register for each (dir, side)
        std::vector<Array<Array<int, 2>, DIM>> m_fineIndices;
        std::vector<Array<Array<Stencil<T>, 2>, DIM>> m_avgFlux;
    }; //end class MBLevelFluxRegister
#include "implem/Proto_MBLevelFluxRegisterImplem.H"
} // end namespace Proto
#endif // end include guard
