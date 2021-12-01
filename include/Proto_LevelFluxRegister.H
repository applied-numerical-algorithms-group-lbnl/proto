#pragma once
#ifndef _PROTO_FLUX_REGISTER_
#define _PROTO_FLUX_REGISTER_

#include "Proto.H"
#include "Proto_MemType.H"
#include "Proto_Register.H"
#include "Proto_BoxData.H"
#include "Proto_Copier.H"
#define LFRT "template<typename T, unsigned int C, MemType MEM>"
#define LFR "LevelFluxRegister<T,C,MEM>"
namespace Proto
{
    // Forward Declarations
    template<typename T, unsigned int C, MemType MEM>
    class LevelFluxRegister;
    template<typename T, unsigned int C, MemType MEM>
    class FluxRegisterCopier;

    // =======================================================================
    // FLUX REGISTER COPIER OP

    /// Flux Register Copier Op
    /**
     */
    template<typename T, unsigned int C=1, MemType MEM = MEMTYPE_DEFAULT> 
    class FluxRegisterCopierOp
    {
        friend class FluxRegisterCopier<T,C,MEM>;
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
    template<typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT>
    class FluxRegisterCopier : public Copier<FluxRegisterCopierOp<T, C, MEM> >
    {
        public:

        inline FluxRegisterCopier(){};

        /// Build Copier Motion Plan
        inline void buildMotionPlans(FluxRegisterCopierOp<T, C, MEM>& a_op);
    }; //end class FluxRegisterCopier

    // =======================================================================
    // FLUX REGISTER

    /// Flux Register
    /**
     */
    template<typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT>
    class LevelFluxRegister
    {
        friend class FluxRegisterCopierOp<T,C,MEM>;
        friend class FluxRegisterCopier<T,C,MEM>;

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
        void incrementCoarse(
                const BoxData<T,C,MEM>& a_flux,
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

        /// Reset
        void reset();

        /// compute the hash key for a coarse register.
        int key(const Box& a_bx,const DataIndex& a_di);

        void print(std::string a_name) const;
        
        protected:
        
        BoxData<T,C,MEM>& sourceData(const Box& a_bx,const DataIndex& di);
        BoxData<T,C,MEM>& destData(const Box& a_bx,const DataIndex& di);

        DisjointBoxLayout  m_crseLayout;
        DisjointBoxLayout  m_fineLayout;
        FluxRegisterCopier<T,C,MEM> m_copier;
        Point              m_refRatio;

        vector<vector<Register<T,C,MEM> > >              m_crseRegisters;
        vector<vector<shared_ptr<BoxData<T,C,MEM> > > >  m_crseIncrement;
        vector<vector<Register<T,C,MEM> > >              m_fineRegisters; 
        // crseIndices is a map because there are in general multiple registers for each (dir, side)
        vector<unordered_map<unsigned int,unsigned int> >m_crseIndices;
        // fineIndices are not a map because there is at most one register for each (dir, side)
        vector<array<array<int,2>, DIM> >                m_fineIndices;
        array<array<Stencil<T>,2>,DIM>                   m_avgFlux;
    }; //end class LevelFluxRegister
#include "implem/Proto_LevelFluxRegisterImplem.H"
} // end namespace Proto
#endif // end include guard
