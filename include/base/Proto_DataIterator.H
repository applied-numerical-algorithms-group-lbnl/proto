#pragma once
#ifndef _PROTO_DATA_ITERATOR_
#define _PROTO_DATA_ITERATOR_
#include <vector>
#include "Proto_DataIndex.H"
//#include "Proto_BoxPartition.H"

namespace Proto
{
    /// Distributed Data Iterator
    /**
        Iterate through the local contents of a DisjointBoxLayout or LevelBoxData.
        Unline normal iterators, DataIterator dereferences into a DataIndex which can
        be used to index into a DisjointBoxLayout, LevelBoxData, or any other "level scoped"
        data holder.

        The intended syntax for for iterating is:
        @code
        DisjointBoxLayout layout(...);
        LevelBoxData<double> data(layout, ...);

        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            Box b = layout[*iter];
            BoxData<double>& patch = data[*iter];
            // etc
        }
        // OR
        for (auto iter = layout.begin(); iter != layout.end(); ++iter)
        {
            Box b = layout[*iter];
            BoxData<double>& patch = data[*iter];
            // etc
        }
        @endcode
    */
    template<typename P>
    class DataIterator
    {
        public:
       
        /// Manual Constructor
        /**
            Construct a DataIterator from an existing DisjointBoxLayout.
            This constructor is mostly used internally. Users can rely on 
            the STL style iteration syntax of container objects
            (e.g. LevelBoxData or DisjointBoxLayout)
        */ 
        inline DataIterator(std::shared_ptr<P> a_partition);
        
        /// Set To Start
        /**
            Modifies the *this and returns it
        */
        inline DataIterator<P>& begin();

        /// Set To End 
        /**
            modifies *this and returns it
        */
        inline DataIterator<P>& end();
        
        /// Continue Query
        /**
            Shortcut to check if *this has reached its end.
        */
        inline bool ok() const;

        /// Increment 
        inline DataIterator<P>& operator++();
       
        /// Equality
        inline bool operator==(const DataIterator<P>& a_rhs) const { return m_current == a_rhs.m_current; }
        
        /// Equality
        inline bool operator==(const DataIndex<P>& a_rhs) const { return m_current == a_rhs; }
        
        /// Inquality
        inline bool operator!=(const DataIterator<P>& a_rhs) const { return m_current != a_rhs.m_current; }
        
        /// Inquality
        inline bool operator!=(const DataIndex<P>& a_rhs) const { return m_current != a_rhs; }

        /// Dereference
        /**
            Unlike normal iterators, dereferencing a DataIterator returns a DataIndex.
            DataIndex is accepted by all "level scope" data holders (e.g. LevelBoxData or DisjointBoxLayout).
        */
        inline const DataIndex<P>& operator*() const;

        inline bool compatible(const P& a_partition) const;

        inline bool compatible(const DataIndex<P>& a_index) const;
        
        inline bool compatible(const DataIterator<P>& a_iter) const;
        
        protected:
        
        std::shared_ptr<P> m_partition;
        DataIndex<P> m_current;
        DataIndex<P> m_end;
    };


#include "implem/Proto_DataIteratorImplem.H"
} // end Proto namespace.
#endif

