#pragma once
#ifndef _PROTO_DATA_INDEX_
#define _PROTO_DATA_INDEX_
#include <vector>
#include <iostream>
#include "Proto.H"
//#include "Proto_BoxPartition.H"
namespace Proto
{
    class DisjointBoxLayout;
    template<typename PP> class DataIterator;

    template<typename P>
    class DataIndex
    {
        friend class DisjointBoxLayout;
        template<typename PP> friend class DataIterator;
        friend std::ostream& operator<<(std::ostream& a_os, const DataIndex& a_index);

        public:

        DataIndex(){};
        inline DataIndex(const shared_ptr<P> a_partition, unsigned int a_index)
        {
            m_partition =  a_partition;
            m_currentInt = a_index;     
        };

        ~DataIndex(){};
        
        inline bool operator==(const DataIndex<P>& a_di) const {return (m_partition == a_di.m_partition) && (m_currentInt == a_di.m_currentInt);};
        inline bool operator!=(const DataIndex<P>& a_di) const {return !(*this==a_di);};

        inline operator int() const { return local(); };
    
        inline int local() const { return m_currentInt - m_partition->procStartIndex(procID()); } 
        inline int global() const { return m_currentInt; } 

        // defined in DataIterator's implem due to circular dependency
        inline bool operator==(const DataIterator<P>& a_iter) const;
        inline bool operator!=(const DataIterator<P>& a_iter) const;

        bool compatible(const DataIndex<P>& a_index) const
        {
            return m_partition->compatible(*a_index.m_partition);
        }
        
        bool compatible(const P& a_partition) const
        {
            return m_partition->compatible(a_partition);
        }

        protected:

        shared_ptr<P> m_partition;
        int m_currentInt;

    }; // end class DataIndex;
    
    template<typename P>
    inline std::ostream& operator<<(std::ostream& os, const DataIndex<P>& a_di)
    {
        os << "DataIndex: partition pointer: " << a_di.m_partition;
        os <<  " | global(local) index: " << a_di.global();
        os << "(" << a_di.local() << ")" << std::endl;
        return os;
    }
}
#endif
