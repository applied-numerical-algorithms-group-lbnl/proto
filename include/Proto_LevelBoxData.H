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
    
    /// Level Copier Operator
    /**
        Copier operator which provides the necessary interface to copy to and from
        LevelBoxData using the Copier interface. See Proto_Copier.H.
    */
    template
    <typename T, unsigned int C=1, MemType MEM = MemType::HOST ,unsigned int CTR=0 >
    class LevelCopierOp
    {
        public:
        LevelBoxData<T,C,MEM,CTR>*  m_src;
        LevelBoxData<T,C,MEM,CTR>*  m_dst;

        inline LevelCopierOp(){};
        
        inline LevelCopierOp(LBD& a_src, LBD& a_dst);
        
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
    };
// =======================================================================
// LEVEL BOX DATA COPIERS
    
    /// Level Copier
    /**
        A Copier used to implement copying of data between LevelBoxData objects.
    */
    template<typename T, unsigned int C=1, MemType MEM=MemType::HOST, unsigned int CTR=0>
    class LevelCopier : public Copier<LevelCopierOp<T, C, MEM, CTR>>
    {
        public:
        
        inline LevelCopier(){};
        
        /// Build Copier Motion Plan
        inline void buildMotionPlans(LevelCopierOp<T, C, MEM, CTR>& a_op);
    };
    
    /// Exchange Copier
    /**
        A Copier used to execute the data copies necessary for filling ghost regions
        within a LevelBoxData.
    */
    template<typename T, unsigned int C=1, MemType MEM=MemType::HOST, unsigned int CTR=0>
    class LevelExchangeCopier : public Copier<LevelCopierOp<T, C, MEM, CTR>>
    {
        public:

        inline LevelExchangeCopier(){};
        
        /// Build Copier Motion Plan
        inline void buildMotionPlans(LevelCopierOp<T, C, MEM, CTR>& a_op);
    };

// =======================================================================
// LEVEL BOX DATA
    
    /// Level Box Data
    /**
        A distributed data holder implementing the abstraction of a union of logically rectangular patches.
    */
    template
    <typename T, unsigned int C=1, MemType MEM=MEMTYPE_DEFAULT, unsigned int CTR=DIM >
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

        /// Function Initialization
        /**
            Initializes *this using a function with a signature similar to that used by
            Proto::forallInPlace_p. This function is purely for convenience. 

            Func signature:
            void a_func(Point& a_point, Var<T, C, MEM>& a_thisData, ... )

            Unlike Proto::forall and its variants, this function cannot manage additional
            data holders (e.g. Srcs should not include a LevelBoxData or similar complex object).
            Elements of a_srcs should be plain-old-data such as single values, Points, or Boxes. 
        */
        template<typename Func, typename... Srcs>
        inline void initialize(Func& a_func, Srcs... a_srcs);

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
        
        /// Linear Offset
        /**
          Returns the serial index of the first data element of this stored on
          this a_proc.
         */
        inline unsigned int offset(int a_proc) const;
        inline unsigned int offset() const { return offset(Proto::procID()); }

        /// Patch Size
        /**
          Computes the size (in data cells) of a single patch in *this.
          The resulting patch ACCOUNTS FOR GHOST CELLS
         */
        inline unsigned int patchSize() const;

        /// Linear In (Patch, Box)
        /**
          Reads data from a serial buffer, populating the patch associated with 
          a_index within the subset specified by a_box. 

          \param a_buffer     The read buffer
          \param a_box        A box defining the subset to read into
          \param a_index      The index defining the patch to read into
         */ 
        inline void linearIn(void* a_buffer, 
                const Box& a_box,
                const DataIndex& a_index);

        /// Linear Out (Patch, Box)
        /**
          Writes data to a serial buffer, from the patch associated with 
          a_index within the subset specified by a_box. 

          \param a_buffer     The write buffer
          \param a_box        A box defining the subset to write from
          \param a_index      The index defining the patch to write from
         */ 
        inline void linearOut(void* a_buffer, 
                const Box& a_box,
                const DataIndex& a_index) const;

        /// Serial Size (Patch, Box)
        /**
          Computes the size in bytes of the serial data buffer needed to hold
          the patch associated with a_index within the subset specified by a_box.

          \param a_box        A box defining the subset of a patch
          \param a_index      The index defining the patch
         */ 
        inline size_t linearSize(const Box& a_box,
                const DataIndex& a_index) const;
        /// Linear In (Patch)
        /**
          Overload of LinearIn that always reads in a full patch
         */
        inline void linearIn(void* a_buffer, 
                const DataIndex& a_index);
        /// Linear Out (Patch)
        /**
          Overload of LinearOut that always writes out a full patch
         */
        inline void linearOut(void* a_buffer, 
                const DataIndex& a_index) const;
        /// Size (Patch)
        /**
          Overload of Size that always computes the size of a full patch
         */
        inline size_t linearSize( const DataIndex& a_index) const;

        /// Linear In (All Local Data)
        /**
          Overload of LinearIn that reads all local data
         */
        inline void linearIn(void* a_buffer);
        /// Linear Out (All Local Data)
        /**
          Overload of LinearOut that writes out all local data
         */
        inline void linearOut(void* a_buffer) const;
        /// Size (All Local Data)
        /**
          Overload of Size that computes the full local size of the data in *this
         */
        inline size_t linearSize() const;

        /// Get Box Layout
        inline DisjointBoxLayout layout() const { return m_dbl; }
        /// Get Box Layout (Overload)
        inline DisjointBoxLayout getDBL() const { return m_dbl; }

        /// Get Ghost Size
        inline Point ghost() const { return m_ghost; }
        /// Get Ghost Size (Overload)
        inline Point getGhostVector() const {return m_ghost;}
        
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
