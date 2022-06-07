#pragma once
#ifndef _PROTO_DATA_INDEX_
#define _PROTO_DATA_INDEX_
#include <vector>
#include <iostream>
#include "Proto.H"
#include "Proto_BoxPartition.H"
namespace Proto
{
    class DisjointBoxLayout;
    class DataIterator;
    class DataIndex
    {
        friend class DisjointBoxLayout;
        friend class DataIterator;
        friend std::ostream& operator<<(std::ostream& a_os, const DataIndex& a_index);

        public:

        DataIndex(){};
        inline DataIndex(const shared_ptr<BoxPartition> a_partition, unsigned int a_index)
        {
            m_partition =  a_partition;
            m_currentInt = a_index;     
        };

        ~DataIndex(){};
        
        inline bool operator==(const DataIndex& a_di) const {return (m_partition == a_di.m_partition) && (m_currentInt == a_di.m_currentInt);};
        inline bool operator!=(const DataIndex& a_di) const {return !(*this==a_di);};

        inline operator int() const { return local(); };
    
        inline int local() const { return m_currentInt - m_partition->procStartIndex(procID()); } 
        inline int global() const { return m_currentInt; } 

        // defined in DataIterator's implem due to circular dependency
        inline bool operator==(const DataIterator& a_iter) const;
        inline bool operator!=(const DataIterator& a_iter) const;

        bool compatible(const DataIndex& a_index) const
        {
            return m_partition->compatible(*a_index.m_partition);
        }
        
        bool compatible(const BoxPartition& a_partition) const
        {
            return m_partition->compatible(a_partition);
        }

        protected:

        shared_ptr<BoxPartition> m_partition;
        int m_currentInt;

    }; // end class DataIndex;
    
    inline std::ostream& operator<<(std::ostream& os, const DataIndex& a_di)
    {
        os << "DataIndex: partition pointer: " << a_di.m_partition;
        os <<  " | global(local) index: " << a_di.global();
        os << "(" << a_di.local() << ")" << std::endl;
        return os;
    }
}
#endif
