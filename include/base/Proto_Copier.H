#pragma once
#ifndef _PROTO_COPIER_
#define _PROTO_COPIER_

#include "Proto_MemType.H"
#include "Proto_Memory.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto.H"
#include "Proto_SPMD.H"
#include "Proto_DataIterator.H"
#include <unordered_map>
#include <cstdint>
#include "Proto_MayDay.H"

namespace Proto

{
// =======================================================================
// DECLARATIONS

    struct BufferEntry;
    struct MotionItem;
    class CopierIterator;
   
// =======================================================================
// ABSTRACT COPIER CLASS

    /// Abstract Generic Parallel Copier
    /**
        Copier is an abstract super-class that handles the MPI communication
        component of a generic parallel copy operation. Copying is executed
        by processing a set of "motion plans" which define how data moves
        from, to, and within this process. 

        The construction of the motion plans must be defined using the abstract
        buildMotionPlans() function. 
    
        The template parameter class OP defines the serialization of any underlying
        dataholders. Generally OP must have a reference to the dataholder(s) themselves
        and must also implement the following API:

        void linearIn (void* a_buffer, Box a_range,  DataIndex a_patchIndex);
        void linearOut(void* a_buffer, Box a_domain, DataIndex a_patchIndex);
        void size(Box a_region, DataIndex a_patchIndex);
        void localCopy(Box a_domain, DataIndex a_domainIndex, Box a_range, DataIndex a_rangeIndex);
    
        linearIn:  Reads data from a_buffer into the part of a patch designated by a_patchIndex intersecting a_range
        linearOut: Writes data to a_buffer from the part of a patch designated by a_patchIndex intersecting a_domain
        size:      Computes the number of elements in the patch designated by a_patchIndex intersectiong a_region
        localCopy: Executes the copy operation between two data patches located on the same parallel process
    */
    
    template<class OP, MemType SRC_MEM, MemType DST_MEM>
    class Copier
    {
        friend class MotionItem;
        friend class CopierIterator; 
        
        public:
        
        enum MotionType
        {
            LOCAL,
            FROM,
            TO
        };

        Copier();
        Copier(OP a_op);
        Copier(Copier& a_copier) = delete;

        // Destructor is virtual to handle potential polymorphic destruct
        inline virtual ~Copier();
        inline void define(OP a_op); 
        inline void clear(); 
        virtual void buildMotionPlans(OP& a_op) = 0;

        inline Copier<OP, SRC_MEM, DST_MEM>& operator=(const Copier<OP, SRC_MEM, DST_MEM>& a_rhs) = delete;
        inline bool operator==(const Copier<OP, SRC_MEM, DST_MEM>& a_rhs) const;
        inline void execute();
        inline void sort();
        inline CopierIterator begin(MotionType a_type) const;
        inline void print() const;
        
        protected:

        OP m_op;

        std::vector<MotionItem> m_localMotionPlan;
        std::vector<MotionItem> m_fromMotionPlan;
        std::vector<MotionItem> m_toMotionPlan;

        bool m_isDefined = false;
        
        private: 
       
        void makeItSo();
        void makeItSoBegin();
        void makeItSoLocal();
        void makeItSoEnd();

        void completePendingSends() const;
        void allocateBuffers();
        void writeToSendBuffers();
        void postSends() const;
        void readFromRecvBuffers();
        void postRecvs() const;

        // Copy buffers        
        void clearBuffers();
        void clearMotionPlan();

        size_t m_sendCapacity;
        size_t m_recvCapacity;
        void* m_sendBuffer;
        void* m_recvBuffer;


#ifdef PR_MPI
        mutable std::vector<BufferEntry> m_fromMe;
        mutable std::vector<BufferEntry> m_toMe;
        
        mutable int m_numSends;
        mutable int m_numRecvs;
        
        mutable std::vector<MPI_Request> m_sendRequests, m_recvRequests;
        mutable std::vector<MPI_Status>  m_sendStatus,  m_recvStatus;
#endif
    }; // end class Copier
    
// =======================================================================
// STRUCTURES

    struct BufferEntry
    {
        void* bufferPointer;
        size_t size;
        const MotionItem* item;
        unsigned int procID;

        inline bool operator<(const BufferEntry& a_rhs) const;
    }; // end struct BufferEntry
    
    struct MotionItem
    {
        DataIndex fromIndex, toIndex;
        Box fromRegion;
        Box toRegion;
        int procID;

        inline MotionItem(){}

        inline MotionItem(const DataIndex& a_fromIndex,
                const DataIndex& a_toIndex,
                const Box& a_fromRegion,
                const Box& a_toRegion);

        inline bool operator==(const MotionItem& a_rhs) const;
        inline bool operator<(const MotionItem& a_rhs) const;
    }; // end struct MotionItem
    
    /*
    struct MotionItemSorter
    {
        inline bool operator()(MotionItem* const & a_lhs,
                MotionItem* const& a_rhs) const
                {
                    return ((*a_lhs) < (*a_rhs));
                }
    }; // end class MotionItemSorter
    */
// =======================================================================
// COPIER ITERATOR

    class CopierIterator
    {
        public:

            inline CopierIterator(const std::vector<MotionItem>* a_motionPlan);
            inline const MotionItem& operator()() const;
            inline const MotionItem& operator[](size_t a_index) const;
            inline const MotionItem& operator*() const;
            inline void operator++();
            inline bool ok(); 
            inline void reset();
            inline size_t size() const;

        private:

            const vector<MotionItem>* m_motionPlan;
            unsigned int m_current;
    }; // end class CopierIterator
    

#include "implem/Proto_CopierImplem.H"
} // end namespace Proto
#endif //end include guard
