#pragma once
#ifndef _PROTO_MMB_PROBLEM_DOMAIN_
#define _PROTO_MMB_PROBLEM_DOMAIN_

#include "Proto_ProbemDomain.H"
#include "Proto_BlockIterator.H"

namespace Proto
{

    /// Mapped Multi-Block Problem Domain
    /**
        The ProblemDomain equivalent for a MMB domain. Contains all of the
        normal ProblemDomain objects for each individual block as well
        as the connectivity information between blocks.
    */
    class MMBProblemDomain
    {
        public:
        
        /// Default Constructor
        inline MMBProblemDomain();
        
        // A NOTE ON CONSTRUCTION
        // The argument list for making an entire MMB ProblemDomain is too deep.
        // Instead, we build up the object incrementally using addBlock and addBoundary.

        /// Add Block
        /**
            Adds a new block to the domain and return the associated index.
            This block is unconnected to other blocks until addBlockBoundary is used
            to add a connection.
        */
        inline BlockIndex addBlock(ProblemDomain& a_domain);

        /// Add Boundary
        /**
            Adds connectivity between two blocks that already exist in *this.
            Returns true if the operation was successful. This function only adds
            boundaries of codimension 1; higher codimensional connectivity are
            added automatically.
            TODO: Need some kind of input that defines the rotation
        */
        inline bool addBlockBoundary(
            const BlockIndex& a_b0,
            const BlockIndex& a_b1);
       
        /// Access Block Problem Domain 
        inline ProblemDomain& operator[](const BlockIndex& a_index);

        // Iterator interface
        inline MMBIterator begin() const;
        inline MMBIterator end() const;

        // Utility
        inline int size() const; 
        inline void print() const;

        private:
        
        // [blockIndex]
        std::vector<ProblemDomain> m_domains; //[blockIndex]
        // Here lies all of the connectivity information
        std::vector<std::array<std::vector<BlockBoundary>, DIM>> m_boundaries; // [blockIndex][codim][boundary]


    }; // end class MMBProblemDomain
}

#include "implem/Proto_MMBProblemDomainImplem.H"
#endif // end include guard
