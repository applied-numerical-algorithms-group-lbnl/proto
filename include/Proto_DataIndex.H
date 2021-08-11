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
    class DataIndex
    {
        protected:
        
        shared_ptr<DBLInternals > m_dblInternals;
        int m_currentInt;
        
        public:
        
        DataIndex(){};
        inline DataIndex(const shared_ptr<DBLInternals> a_dblInternals, int a_index)
        {
            m_dblInternals =  a_dblInternals;
            m_currentInt = a_index;     
        };
        inline void defineEndIndex(const shared_ptr<DBLInternals> a_dblInternals, int a_end)
        {
            m_dblInternals =  a_dblInternals;
            m_currentInt = a_end;
        }
        ~DataIndex(){};
        inline int intIndex() const { return m_currentInt; };
        
        inline bool operator==(const DataIndex& a_di) const
        {return (m_dblInternals == a_di.m_dblInternals) && (m_currentInt == a_di.m_currentInt);};
        
        inline bool operator!=(const DataIndex& a_di) const
        {return !(*this==a_di);};

        inline operator int() const { return m_currentInt; };
    };
    inline std::ostream& operator<<(std::ostream& os, const DataIndex& a_di)
    {
        return os <<  "index " << a_di.intIndex() << " " ;
    }
}
#endif
