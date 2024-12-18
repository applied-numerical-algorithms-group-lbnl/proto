#pragma once
#ifndef _PROTO_DISTRIBUTED_SET_
#define _PROTO_DISTRIBUTED_SET_

#include "Proto.H"

namespace Proto
{
    template<typename DATA>
    class DistributedSet
    {
        public:
        DistributedSet();
        ~DistributedSet();
        void insert(const DATA& item);
        void exchange();
        void clearLocal();
        void clearGlobal();
        std::set<DATA> data() { return m_globalSet; }

        private:
        
        std::set<DATA> m_globalSet;
        #ifdef PR_MPI
        std::vector<DATA> m_localBuffer;

        void* m_sendBuffer;
        void* m_recvBuffer;
        size_t m_sendBufferSizeInBytes;
        size_t m_recvBufferSizeInBytes;

        int* m_localBufferSizesInBytes;
        int* m_bufferOffsetsInBytes;

        #endif

    };

    #include "implem/Proto_DistributedSetImplem.H"



}

#endif