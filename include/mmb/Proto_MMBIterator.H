#pragma once
#ifndef _PROTO_MMB_ITERATOR_
#define _PROTO_MMB_ITERATOR_
#include "Proto_ProbemDomain.H"

namespace Proto
{
    
    // Forward declaration
    class MMBProblemDomain;

    /// Multi-Block Index
    /**
        Identifies a block uniquely
        TODO: Add additional functionality here if needed.
    */
    struct BlockIndex
    {
        int index;
    }; // end struct MultiBlockIndex
    
    
    class MMBIterator
    {
    public:

    MMBIterator(MMBProblemDomain a_domain);
    inline const MMBIndex operator*() const;
    private:

    MultiBlockIndex m_index;

    }; // end class MMBIterator
}

#include "implem/Proto_MMBIteratorImplem.H"
#endif // end include guard
