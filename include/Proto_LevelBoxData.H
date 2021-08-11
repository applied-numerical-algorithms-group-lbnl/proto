#pragma once
#ifndef _PROTO_LEVELBOXDATA_H_
#define _PROTO_LEVELBOXDATA_H_

#include "Proto_MemType.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_Copier.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_LevelIterators.H"
#include "Proto_WriteBoxData.H"
#include "Proto_MayDay.H"
#include "Proto_SPMD.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cstring> // for Writing data to .vtk files.
#define LBD LevelBoxData<T, C, MEM, CTR>

//using namespace std;
//using std::shared_ptr;

namespace Proto 
{
    // Forward Declarations
    template<typename T, unsigned int C, MemType MEM, unsigned int CTR>
    class LevelBoxData;

// =======================================================================
// LEVEL BOX DATA COPIER OP
    template
    <typename T, unsigned int C=1, MemType MEM = MemType::HOST ,unsigned int CTR=0 >
    class LevelCopierOp
    {
        public:
        LevelBoxData<T,C,MEM,CTR>*  m_src;
        LevelBoxData<T,C,MEM,CTR>*  m_dst;

        inline LevelCopierOp(){};
        
        inline LevelCopierOp(LBD& a_src, LBD& a_dst);

        inline int linearSize(const Box& a_bx,
                const DataIndex& a_index) const;
        
        inline void linearOut(void* a_buf, 
                const Box& a_bx,
                const DataIndex& a_index) const;
        
        inline void linearIn(void* a_buf, 
                const Box& a_bx,
                const DataIndex& a_index) const;
        
        inline void localCopy(
            const Box& a_domain, const DataIndex& a_domainIndex,
            const Box& a_range,  const DataIndex& a_rangeIndex) const;
    };
// =======================================================================
// LEVEL BOX DATA COPIERS
    
    template<typename T, unsigned int C=1, MemType MEM=MemType::HOST, unsigned int CTR=0>
    class LevelCopier : public Copier<LevelCopierOp<T, C, MEM, CTR>>
    {
        public:
        
        inline LevelCopier(){};
        inline void buildMotionPlans(LevelCopierOp<T, C, MEM, CTR>& a_op);
    };
    
    template<typename T, unsigned int C=1, MemType MEM=MemType::HOST, unsigned int CTR=0>
    class LevelExchangeCopier : public Copier<LevelCopierOp<T, C, MEM, CTR>>
    {
        public:

        inline LevelExchangeCopier(){};
        inline void buildMotionPlans(LevelCopierOp<T, C, MEM, CTR>& a_op);
    };

// =======================================================================
// LEVEL BOX DATA
    
    /// Level Box Data
    /**
      A distributed data holder implementing the abstraction of a union of logically rectangular patches.
     */
    template
    <typename T, unsigned int C=1, MemType MEM=MemType::HOST, unsigned int CTR=DIM >
    class LevelBoxData
    {
        public:

        /// Default constructor
        /**
            Constructs an undefined LevelBoxData. Can be initialized through define(...)
        */
        inline LevelBoxData();

        /// Layout Constructor
        /**
            Constructs a non-trivial LevelBoxData by calling define.
        */
        inline LevelBoxData(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

        /// Define
        /**
            Lazily constructs a non-trivial LevelBoxData
        */
        inline void define(const DisjointBoxLayout& a_dbl, const Point& a_ghost);

        /// Size
        /**
            Returns the number of boxes in the union across all processes.
        */
        inline unsigned int size() const;

        /// Initialize To Zero
        /**
            Sets all values in all patches equal to zero. This function is mostly for 
            debugging purposes and should be used sparingly if at all in production code.
        */
        inline void setToZero();

        /// Exchange
        /**
            Copies data from the valid regions of *this into ghost regions. When MPI is
            enabled, this function also takes care of any necessary communication between
            patches on different processes.
        */
        inline void exchange();

        /// Copy To
        /**
            Copies data from the valid regions of this into the valid regions of another
            LevelBoxData. This function only operates on VALID data. To fill ghost regions,
            use exchange().
        */
        inline void copyTo(LevelBoxData& a_dest) const;

        /// Get Box Layout
        inline DisjointBoxLayout layout(){ return m_dbl; }
        /// Get Box Layout (Overload)
        inline DisjointBoxLayout getDBL(){ return m_dbl; }

        /// Get Ghost Size
        inline Point ghost() { return m_ghost; }
        /// Get Ghost Size (Overload)
        inline Point getGhostVector() {return m_ghost;}
        
        //=======================================
        // DATA ITERATOR INTERFACE

        /// Get Patch
        inline BoxData<T,C,MEM>& operator[](const DataIndex& a_index);
        /// Get Const Patch
        inline const BoxData<T,C,MEM>& operator[](const DataIndex& a_index) const;
        
        /// Get Iterator
        inline DataIterator begin() const;

        /// For debugging purposes.
        int s_verbosity = 0;
        
        private: 

        std::vector<shared_ptr<BoxData<T,C> > > m_data;
        Point                                   m_ghost;
        DisjointBoxLayout                       m_dbl;
        bool                                    m_isDefined;
        LevelExchangeCopier<T, C, MEM, CTR>     m_exchangeCopier;
        
        // disallow copy constructors and assignment operators
        // to avoid very hard-to-find performance problems
        LevelBoxData<T, C, MEM, CTR>& operator= (const LevelBoxData<T, C, MEM, CTR>& a_rhs);
        LevelBoxData<T, C, MEM, CTR>(const LevelBoxData<T, C, MEM, CTR>& a_rhs);
    }; // end class LevelBoxData
    

#include "implem/Proto_LevelBoxDataImplem.H"
} //end namespace Proto
#endif
