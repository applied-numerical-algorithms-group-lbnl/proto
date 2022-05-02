#pragma once
#ifndef _PROTO_BLOCK_ITERATOR_
#define _PROTO_BLOCK_ITERATOR_
#include "Proto_ProbemDomain.H"
#include "Proto_BlockIndex.H"

namespace Proto
{
    
    // Forward declaration
    class MMBProblemDomain;
    class MMBIterator
    {
    public:

    BlockIterator(MMBProblemDomain a_domain);
    inline const BlockIndex operator*() const;
    private:

    BlockIndex m_index;

    }; // end class MMBIterator
}

#include "implem/Proto_MMBIteratorImplem.H"
#endif // end include guard
