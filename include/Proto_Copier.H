#pragma once

#include "Proto_DisjointBoxLayout.H"
#include "Proto.H"
#include "Proto_SPMD.H"
#include "Proto_LevelIterators.H"
#include <unordered_map>
#include <cstdint>

namespace Proto
{
    template<class OP>
    class Copier
    {
        friend class MotionItem;
        friend class CopyIterator; 
        public:

        Copier():m_isDefined(false){}
        Copier(const Copier<OP>& a_rhs);

        // Destructor is virtual to handle potential polymorphic destruct
        inline virtual ~Copier();
        
        /// Define
        /**
            This function MUST be implemented by classes that inherit
            from Copier in order to specify the construction of the
            motion plans.
        */
        inline virtual void define() = 0;
        inline virtual void clear(); 

        inline Copier<OP>& operator=(const Copier<OP>& a_rhs);
        inline bool operator==(const Copier<OP>& a_rhs) const;
        inline void reverse();
        inline void execute();
        
        protected:
        

        OP* m_op;

        std::vector<MotionItem*> m_localMotionPlan;
        std::vector<MotionItem*> m_fromMotionPlan;
        std::vector<MotionItem*> m_toMotionPlan;

        bool m_isDefined;
        inline void sort();
        
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
        
    }; // end class Copier
    
    struct BufferEntry
    {
        void* bufferPointer;
        size_t size;
        const MotionItem* item;
        unsigned int procID;

        bool operator<(const BufferEntry& a_rhs) const;
    }; // end struct BufferEntry
    
    class MotionItem
    {
        public:

            DataIndex fromIndex, toIndex;
            Box fromRegion;
            Box toRegion;
            int procID;

            MotionItem(const DataIndex& a_fromIndex,
                    const DataIndex& a_toIndex,
                    const Box& a_fromRegion,
                    const Box& a_toRegion);

            bool operator==(const MotionItem& a_rhs) const;
    }; // end class MotionItem
    
    class CopyIterator
    {
        public:
            enum MotionType
            {
                LOCAL,
                FROM,
                TO
            }

            inline CopyIterator(const Copier& a_copier, MotionType);
            inline const MotionItem& operator()() const;
            inline const MotionItem& operator[](size_t a_index) const;
            inline void operator++();
            inline void ok(); 
            inline void reset();
            inline size_t size();

        private:

            const vector<MotionItem*>* m_motionPlan;
            unsigned int m_current;
    }; // end class CopyIterator
#include "implem/Proto_CopierImplem.H"
} // end namespace Proto
