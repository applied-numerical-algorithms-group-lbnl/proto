#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef _Proto_DataIndex_
#define _Proto_DataIndex_
#include <vector>
#include "Proto.H"
#include "Proto_DBLInternals.H"
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
        inline DataIndex(const shared_ptr<DBLInternals> a_internals, int a_index)
        {
            m_internals =  a_internals;
            m_currentInt = a_index;     
        };

        ~DataIndex(){};
        
        inline bool operator==(const DataIndex& a_di) const {return (m_internals == a_di.m_internals) && (m_currentInt == a_di.m_currentInt);};
        inline bool operator!=(const DataIndex& a_di) const {return !(*this==a_di);};

        inline operator int() const { return local(); };
    
        inline int local() const { return m_currentInt - m_internals->startProc[Proto::procID()]; } 
        inline int global() const { return m_currentInt; } 

        // defined in DataIterator's implem due to circular dependency
        inline bool operator==(const DataIterator& a_iter) const;
        inline bool operator!=(const DataIterator& a_iter) const;

        bool compatible(const DataIndex& a_index) const
        {
            return m_internals->compatible(*a_index.m_internals);
        }
        
        bool compatible(const DBLInternals& a_internals) const
        {
            return m_internals->compatible(a_internals);
        }

        protected:

        shared_ptr<DBLInternals > m_internals;
        int m_currentInt;

    }; // end class DataIndex;
    
    inline std::ostream& operator<<(std::ostream& os, const DataIndex& a_di)
    {
        os << "DataIndex: internals pointer: " << a_di.m_internals;
        os <<  " | global(local) index: " << a_di.global();
        os << "(" << a_di.local() << ")" << std::endl;
        return os;
    }
}
#endif
