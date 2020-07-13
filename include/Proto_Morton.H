#ifndef _PROTO_MORTON_H_
#define _PROTO_MORTON_H_
#include <array>
#include <cstdint>
#include "Proto_Point.H"
using namespace std;
#if DIM==2
#define LOGSIZE 32 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 4294967296 // = 2^LOGSIZE.
#endif
#if DIM==3
#define LOGSIZE 21 // DIM*LOGSIZE must be less than the size (in bits) of (uint64_t).
#define MORTONSIZE 2097152 // = 2^LOGSIZE.
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
  class MortonOrder
  {
  private:
#if DIM > 1 
#if DIM < 7
    static array<array<uint64_t,MORTONSIZE>,DIM> s_morton1D;
#endif
#endif
    static bool s_isDefined;

    // Static define must be called before using operator(). 
    inline static void define()
    {
      PR_assert(DIM < 7);
      s_isDefined = true;
#if DIM > 1
      uint64_t mask0 = 1;
      for (uint64_t i = 0; i <MORTONSIZE; i++)
        {
          for (uint64_t d = 0; d < DIM;d++)
            {
              s_morton1D[d][i] = 0;
            }
          for (uint64_t logi = 0; logi < LOGSIZE; logi++)
            {
              for (uint64_t d = 0;d < DIM; d++) 
                {
                  s_morton1D[d][i] += 
                    ((i >> logi) & mask0) << (DIM*logi + d);
                }
            }
        }
#endif
    };
    // operator that computes the Morton index corresponding to the input Point.
  public:
    inline static uint64_t index(const Point& a_pt)
    {
      if (!s_isDefined) define();
      uint64_t retval = 0;
      // PC::What's wrong ?? PR_assert(a_pt%(MORTONSIZE*Point::Ones()) == a_pt);
      // We are enforcing the condition that the Morton index can fit in a 64-bit integer.
#if DIM > 1
      for (int d = 0; d < DIM; d++)
        {
          retval += s_morton1D[d][a_pt[d]];
          // retval += s_morton1D[d][a_pt[d]%MORTONSIZE]; Used for many-to-one Morton index.
        }
#else
      retval = a_pt[0];
#endif
      return retval;
    };
    inline uint64_t operator()(const Point& a_pt) const {return index(a_pt);}
  };
  class MortonSort
  {
  public:
    MortonSort();
    ~MortonSort();
    inline static bool compare(const Point& pt1,const Point& pt2)
    {
      return (MortonOrder::index(pt1) < MortonOrder::index(pt2));
    };
    inline void sort(vector<Point>& a_pts)
    {
      std::sort(a_pts.begin(),a_pts.end(),compare);
    };
  };
} // end Proto namespace.
#endif
