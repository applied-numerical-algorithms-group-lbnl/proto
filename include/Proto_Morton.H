#ifndef _PROTO_MORTON_H_
#define _PROTO_MORTON_H_
#include <array>
#include <cstdint>
#include <algorithm>
#include "Proto_Point.H"
using namespace std;
#if DIM==2
#define LOGSIZE 16 // 21 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 65536 // 2097152 // = 2^LOGSIZE.
#endif
#if DIM==3
#define LOGSIZE 10 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 1024 // = 2^LOGSIZE.
#endif
#if DIM==4
#define LOGSIZE 16 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 65536 // = 2^LOGSIZE.
#endif
#if DIM==5
#define LOGSIZE 12 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 4096 // = 2^LOGSIZE.
#endif
#if DIM==6
#define LOGSIZE 10 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 1024 // = 2^LOGSIZE.
#endif
namespace Proto
{
    /// Morton Indexer
    /** 
       Class for computing Morton Index of a DIM-tuple corresponding to the 
       bits of each element of the tuple.
    */
    class Morton
    {
        public:
#if DIM > 1 
#if DIM < 7
            array<vector<uint64_t>,DIM> m_morton1D;
#endif
#endif
            /// Constructor
            Morton()
            {
                PR_TIMERS("Morton define");
                PR_assert(DIM < 7);
                for (int d = 0; d < DIM; d++)
                {
                    m_morton1D[d]=vector<uint64_t>(MORTONSIZE);
                }
#if DIM > 1
                uint64_t mask0 = 1;
                for (uint64_t i = 0; i <MORTONSIZE; i++)
                {
                    for (uint64_t d = 0; d < DIM;d++)
                    {
                        m_morton1D[d][i] = 0;
                    }
                    for (uint64_t logi = 0; logi < LOGSIZE; logi++)
                    {
                        for (uint64_t d = 0;d < DIM; d++) 
                        {
                            m_morton1D[d][i] += 
                                ((i >> logi) & mask0) << (DIM*logi + d);
                        }
                    }
                }
            };
#endif
    };
    
    /// Compute Morton Index
    /**
        Compute the Morton index of a Point.

        \param a_pt     A Point
        \param a_morton A Morton indexer
    */
    inline uint64_t mortonIndex(const Point& a_pt,const Morton& a_morton)
    {
        uint64_t retval = 0;
#if DIM > 1
        for (int d = 0; d < DIM; d++)
        {
            retval +=m_morton.a_morton1D[d][a_pt[d]];
        }
#else
        retval = a_pt[0];
#endif
        return retval;
    }
    
    // Sorting condition
    inline bool compareSecond(pair<Point,uint64_t> a_a,pair<Point,uint64_t> a_b)
    {
        return a_a.second < a_b.second;
    }  
    
    /// Morton Sort
    /**
        Sort a vector of Points by Morton index.

        \param a_pts    A vector of Points to sort
        \param a_morton A Morton indexer
    */
    inline void mortonSort(vector<Point>& a_pts, const Morton& a_morton)
    {
        vector<pair<Point,uint64_t> > sorter;
        for (int k = 0; k < a_pts.size();++k)
        {
            sorter.push_back(pair<Point,uint64_t>(a_pts[k],mortonIndex(a_pts[k],a_morton)));
        }

        std::sort(sorter.begin(),sorter.end(),compareSecond);
        
        for (int k = 0; k < sorter.size();++k)
        {
            a_pts[k] = sorter[k].first;
        }
    }
}// end Proto namespace.
#endif
