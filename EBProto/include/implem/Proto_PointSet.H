#ifndef _Proto_STDSETIVS_H_
#define _Proto_STDSETIVS_H_

#include "Proto.H"
#include <set>
using std::set;
using std::vector;

namespace Proto
{
  class PointSetIterator;


  /// PointSet implementation based std::set
  /**
   */
  class PointSet
  {
  public:
    friend class PointSetIterator; 

    ///
    PointSet(){;}

    ///
    inline PointSet(const Box& a_box);

    ///
    inline PointSet(const PointSet& a_sivs);

    ///
    ~PointSet(){;}

    ///
    inline void define(const Box&);

    ///
    inline void define(const PointSet& a_sivs);

    ///
    inline PointSet& operator=(const PointSet& a_sivs);

    ///or
    inline PointSet& operator|=(const PointSet& a_sivs);

    ///
    inline PointSet& operator|=(const Point& a_iv);

    ///
    inline PointSet& operator|=(const Box& a_box);

    ///and
    inline PointSet& operator&=(const PointSet& s_sivs);

    ///and
    inline PointSet& operator&=(const Box& a_box);


    ///not
    inline PointSet& operator-=(const PointSet& a_sivs);

    ///not
    inline PointSet& operator-=(const Point& a_iv);

    ///not
    inline PointSet& operator-=(const Box& a_box);

    ///returns true if
    inline bool operator==(const PointSet& lhs) const;

    ///
    inline bool contains(const Point& iv) const;

    ///
    inline bool contains(const Box& box) const;

    ///
    inline void grow(int igrow);

    ///
    inline void grow(int idir, int igrow);

    ///
    inline void refine(int iref = 2);

    ///
    inline void coarsen(int iref = 2);

    ///
    inline void shift(const Point& iv);

    ///
    inline void clear();

    ///
    inline Box minBox() const;

    ///
    inline bool isEmpty() const;

    ///
    inline void makeEmpty() ;

    ///
    inline int numPts() const;

  private:

    std::set<Point> m_stdSet;
  };

  ///
  class PointSetIterator
  {
  public:
    ///
    inline PointSetIterator();

    ///
    inline PointSetIterator(const PointSet& ivs);

    ///
    inline void define(const PointSet& ivs);

    ///
    inline const Point& operator()() const ;

    ///
    inline bool ok() const;

    ///
    inline void operator++();

    ///
    inline void begin();

    ///
    inline void end();

    ///
    inline void clear();
  private:
    const PointSet* m_ivs;
    std::set<Point>::iterator m_iter;
  };
}
#include "Proto_PointSetImplem.H"
#endif 
