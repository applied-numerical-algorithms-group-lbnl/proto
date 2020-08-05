#ifndef _PROTO_MORTON_H_
#define _PROTO_MORTON_H_
#include <array>
#include <cstdint>
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
/* 
Class for computing Morton Index of a DIM-tuple corresponding to the 
bits of each element of the tuple.
*/
namespace Proto
{
  class Morton
  {
  public:
#if DIM > 1 
#if DIM < 7
    array<vector<uint64_t>,DIM> m_morton1D;
#endif
#endif
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
// Function that computes the Morton index corresponding to the input Point.

  inline uint64_t mortonIndex(const Point& a_pt,const Morton& m_morton)
  {
    uint64_t retval = 0;
    PR_assert((a_pt%(Point::Ones()*MORTONSIZE)) == a_pt);
    // We are enforcing the condition that the Morton index can be computed.
#if DIM > 1
    for (int d = 0; d < DIM; d++)
      {
        retval +=m_morton.m_morton1D[d][a_pt[d]];
        // retval += s_morton1D[d][a_pt[d]%MORTONSIZE]; Used for many-to-one Morton index.
      }
#else
    retval = a_pt[0];
#endif
    return retval;
  };
  //static uint64_t operator()(const Point& a_pt) const {return index(a_pt);}
  inline bool compareSecond(pair<Point,uint64_t> a_a,pair<Point,uint64_t> a_b)
    {return a_a.second < a_b.second;};  
  inline void mortonSort(vector<Point>& a_pts, const Morton& m_morton)
  {
    vector<pair<Point,uint64_t> > sorter;
    for (int k = 0; k < a_pts.size();++k)
      {
        sorter.push_back(pair<Point,uint64_t>
                         (a_pts[k],mortonIndex(a_pts[k],m_morton)));
      }
       
    std::sort(sorter.begin(),sorter.end(),compareSecond);
    for (int k = 0; k < sorter.size();++k)
      {
        a_pts[k] = sorter[k].first;
      }
  };
}// end Proto namespace.
#endif
