#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef Proto_LEVELITERATORS_H_
#define Proto_LEVELITERATORS_H_
#include <vector>
#include "Proto.H"
#include "Proto_DataIndex.H"
#include "Proto_DisjointBoxLayout.H"

namespace Proto
{
    /// Distributed Data Iterator
    /**
        Iterate through the local contents of a DisjointBoxLayout or LevelBoxData.
        Unline normal iterators, DataIterator dereferences into a DataIndex which can
        be used to index into a DisjointBoxLayout, LevelBoxData, or any other "level scoped"
        data holder.

        The intended syntax for for iterating is:

        DisjointBoxLayout layout(...);
        LevelBoxData<double> data(layout, ...);

        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            Box b = layout[*iter];
            BoxData<double>& patch = data[*iter];
            // etc
        }
    */
    class DataIterator
    {
        protected:
        
        DisjointBoxLayout m_layout;
        DataIndex m_current;
        DataIndex m_end;
        
        public:
       
        /// Manual Constructor
        /**
            Construct a DataIterator from an existing DisjointBoxLayout.
            Generally there is no reason to do this; it's usually easier to call
            DisjointBoxLayout::begin() or LevelBoxData::begin(). 
            FIXME: Currently the iterator that is constructed is initialized as end() and I'm not sure why... -CLG
        */ 
        inline DataIterator(const DisjointBoxLayout& a_dbl);
        
        /// Set To Start
        /**
            Modifies the *this and returns it
        */
        inline DataIterator& begin();

        /// Set To End 
        /**
            modifies *this and returns it
        */
        inline DataIndex& end(){return m_end;};
        
        /// Continue Query
        /**
            Shortcut to check if *this has reached its end.
        */
        inline bool ok() const;

        /// Increment 
        inline DataIterator& operator++();
        
        /// Dereference
        /**
            Unlike normal iterators, dereferencing a DataIterator returns a DataIndex.
            DataIndex is accepted by all "level scope" data holders (e.g. LevelBoxData or DisjointBoxLayout).
        */
        inline const DataIndex& operator*() const;

        /// Get Patch Point
        /**
            Returns a Point associated with this patch in the level.
        */
        inline Point point() const;

        /// Returns the Box corresponding to that index. Undefined if invalid.
        /**
            Returns the Box associated with this patch in the layout.
            Does not take any possible LevelBoxData ghost cells into account. 
        */
        inline Box box() const;

        inline const DisjointBoxLayout& layout() const;

        inline bool compatible(const DataIndex& a_index) const;
        
        inline bool compatible(const DataIterator& a_iter) const;
    };


    class NeighborIterator
    {
        protected:
        
        DisjointBoxLayout m_layout;
        Box m_bx;
        BoxIterator m_bxit;
        Point m_shift;
        DataIndex m_current;
        DataIndex m_end;
        ProblemDomain m_dProbC;
        
        public:
        
        /// Constructor takes DBL, and Box that we are intersecting with the DBL.
        inline NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx);
        //~NeighborIterator();

        /// Iterator functions begin, end, ++.
        inline NeighborIterator& begin();
        inline const DataIndex& end() const{return m_end;};
        inline bool ok() const {return (m_current != m_end);}
        inline NeighborIterator& operator++();

        /// Dereference operator returns the DataIndex corresponding to the current iterate.
        inline const DataIndex& operator*() const{return m_current;};
        /// The Point corresponding to the Box for that index. Undefined if invalid.
        inline Point point() const;
        /// The region to be copied from corresponding to that index. Undefined if invalid.
        inline Box srcBox() const;
        /// The region to be copied from corresponding to that index. Undefined if invalid.
        inline Box destBox() const;
        /// The amount by which the srcBox is shifted to obtain destBox a_bx.
        inline Point shift();
        // Processor ID for current iterate.
        inline int procID(){return m_layout.procID(m_current);}; 
    };
#include "implem/Proto_LevelIteratorsImplem.H"
} // end Proto namespace.
#endif

