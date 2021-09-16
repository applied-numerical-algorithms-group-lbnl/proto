#include "Proto.H"
#include "Proto_Register.H"
#include "Proto_BoxData.H"
#include "Proto_Copier.H"

namespace Proto
{
    // Forward Declarations
    template<typename T, unsigned int C, MemType MEM>
    class LevelFluxRegister;
    
// =======================================================================
// FLUX REGISTER COPIER OP
    
    /// Flux Register Copier Op
    /**
    */
    template<typename T, unsigned int C=1, MemType MEM = MemType::DEFAULT_MEMTYPE> 
    class FluxRegisterCopierOp
    {
        private:
        LevelFluxRegister<T, C, MEM>* m_register;

        public:
        inline FluxRegisterCopierOp(){};
        inline FluxRegisterCopierOp(LevelFluxRegister<T, C, MEM>& a_register);
        
        /// Get Serialized Size
        inline int linearSize(const Box& a_bx,
                const DataIndex& a_index) const;
        
        /// Serial Write To Buffer
        inline void linearOut(void* a_buf, 
                const Box& a_bx,
                const DataIndex& a_index) const;
        
        /// Serial Read From Buffer
        inline void linearIn(void* a_buf, 
                const Box& a_bx,
                const DataIndex& a_index) const;
        
        /// Local Copy
        inline void localCopy(
            const Box& a_domain, const DataIndex& a_domainIndex,
            const Box& a_range,  const DataIndex& a_rangeIndex) const;

    }; // end class FluxRegisterCopierOp
    
// =======================================================================
// FLUX REGISTER COPIER
    
    /// Flux Register Copier
    /**
    */
    template<typename T, unsigned int C=1, MemType MEM=MemType::DEFAULT_MEMTYPE>
    class FluxRegisterCopier : public Copier<FluxRegisterCopierOp<T, C, MEM>>
    {
        public:
        
        inline FluxRegisterCopier(){};
        
        /// Build Copier Motion Plan
        inline void buildMotionPlans(FluxRegisterCopierOp<T, C, MEM, CTR>& a_op);
    }; //end class FluxRegisterCopier

// =======================================================================
// FLUX REGISTER

    /// Flux Register
    /**
    */
    template<typename T, unsigned int C, MemType MEM>
    class LevelFluxRegister
    {
        friend class FluxRegisterCopierOp;

        public:

        /// Constructor
        LevelFluxRegister(const DisjointBoxLayout& a_dblCoarse,
                const DisjointBoxLayout& a_dblFine,
                const Point& a_refRatio);

        /// Define function
        void define(const DisjointBoxLayout& a_dblCoarse,
                const DisjointBoxLayout& a_dblFine,
                const Point& a_refRatio);

        /// Increment Coarse
        /**
            Increments all the registers contained in the coarse patch a_di, direction a_dir. 
            Assumes that a_flux is defined on all the a_dir faces in the patch.
         */
        void incrementCoarse(const BoxData<T,C,MEM>& a_flux,
                const DataIndex& a_di,
                const T& a_weight,
                unsigned int a_dir);

        /// Increment Fine
        /**
            Increments all the registers contained in the fine patch a_di, direction a_dir. 
            Assumes that a_flux is defined on all the a_dir faces in the patch.
        */
        void incrementFine(const BoxData<T,C,MEM>& a_flux,
                const DataIndex& a_di,
                const T& a_weight,
                unsigned int a_dir); 

        /// Reflux
        /**
            Performs the reflux correction on input coarse data. 
            FIXME: document the exact usage of weight
        */
        void reflux(LevelBoxData<T,C,MEM>& a_coarseData,
                const T& a_weight);

        /// access to fine register data. Used in LFROperator for Copier.
        BoxData<T,C,MEM>& sourceBD(const Box& a_bx,const DataIndex& di);
        
        /// Access to coarse register data. Used in LFROperator for Copier.
        BoxData<T,C,MEM>& destBD(const Box& a_bx,const DataIndex& di);
        
        /// Reset
        void reset();
        
        /// compute the hash key for a coarse register.
        int key(const Box& a_bx,const DataIndex& a_di);

        protected:
        
        DisjointBoxLayout  m_crseLayout;
        DisjointBoxLayout  m_fineLayout;
        FluxRegisterCopier m_copier;
        Point              m_refRatio;

        vector<vector<Register<T,C,MEM>>>                m_coarseRegisters;
        vector<vector<Register<T,C,MEM>>>                m_fineRegisters; 
        vector<vector<shared_ptr<BoxData<T,C,MEM>>>>     m_destData;
        vector<unordered_map<unsigned int,unsigned int>> m_coarseIndices;
        array<array<Stencil<T>,2>,DIM>                   m_avgFlux;
        vector<array<array<int,2>, DIM>>                 m_hasFineRegister;
    }; //end class LevelFluxRegister
#include "implem/Proto_LevelFluxRegisterImplem.H"
} // end namespace Proto
