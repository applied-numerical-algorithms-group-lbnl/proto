#pragma once

#include "Proto_DisjointBoxLayout.H"
#include "Proto.H"
#include "Proto_SPMD.H"
#include "Proto_LevelIterators.H"
#include <unordered_map>
#include <cstdint>

namespace Proto

{
    // forward declarations
    struct BufferEntry;
    struct MotionItem;
    class CopierIterator;
   
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
    
    template<class OP>
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
        }

        Copier();
        Copier(OP& a_op);

        // Destructor is virtual to handle potential polymorphic destruct
        inline virtual ~Copier();
        inline void define(OP& a_op);
        inline void clear(); 
        inline virtual void buildMotionPlans(OP& a_op) = 0;

        inline Copier<OP>& operator=(const Copier<OP>& a_rhs);
        inline bool operator==(const Copier<OP>& a_rhs) const;
        inline void reverse();
        inline void execute();
        inline CopierIterator begin(MotionType a_type) const;
        protected:

        OP* m_op;

        std::vector<MotionItem*> m_localMotionPlan;
        std::vector<MotionItem*> m_fromMotionPlan;
        std::vector<MotionItem*> m_toMotionPlan;

        bool m_isDefined;
        
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
        
        mutable void* m_sendBuffer;
        mutable size_t m_sendCapacity;
        mutable void* m_recvBuffer;
        mutable size_t m_recvCapacity;
        
        mutable std::vector<BufferEntry> m_fromMe;
        mutable std::vector<BufferEntry> m_toMe;
        mutable std::vector<std::vector<BufferEntry>> m_toMeUnpack;        
        
        mutable int m_numSends;
        mutable int m_numRecvs;

#ifdef PR_MPI
        mutable std::vector<MPI_Request> m_sendRequest, m_recvRequests;
        mutable std::vector<MPI_Status>  m_sendStatus,  m_recvStatus;
#endif
    }; // end class Copier
    
    struct BufferEntry
    {
        void* bufferPointer;
        size_t size;
        const MotionItem* item;
        unsigned int procID;

        bool operator<(const BufferEntry& a_rhs) const;
    }; // end struct BufferEntry
    
    struct MotionItem
    {
        DataIndex fromIndex, toIndex;
        Box fromRegion;
        Box toRegion;
        int procID;

        MotionItem(const DataIndex& a_fromIndex,
                const DataIndex& a_toIndex,
                const Box& a_fromRegion,
                const Box& a_toRegion);

        bool operator==(const MotionItem& a_rhs) const;
        bool operator<(const MotionItem& a_rhs) const;
    }; // end class MotionItem
    
    class CopierIterator
    {
        public:

            inline CopierIterator(const std::vector<MotionItem*>* a_motionPlan);
            inline const MotionItem& operator()() const;
            inline const MotionItem& operator[](size_t a_index) const;
            inline void operator++();
            inline void ok(); 
            inline void reset();
            inline size_t size();

        private:

            const vector<MotionItem*>* m_motionPlan;
            unsigned int m_current;
    }; // end class CopierIterator
    
#include "implem/Proto_CopierImplem.H"
} // end namespace Proto
