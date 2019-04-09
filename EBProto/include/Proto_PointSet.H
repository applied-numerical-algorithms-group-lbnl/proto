#ifndef _Proto_STDSETIVS_H_
#define _Proto_STDSETIVS_H_

#include "Proto.H"
#include <set>
using std::set;
using std::vector;
using std::unordered_map;
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
    PointSet(const Box& a_box);

    ///
    PointSet(const PointSet& a_sivs);

    ///
    ~PointSet(){;}

    ///
    void define(const Box&);

    ///
    void define(const PointSet& a_sivs);

    ///
    PointSet& operator=(const PointSet& a_sivs);

    ///or
    PointSet& operator|=(const PointSet& a_sivs);

    ///
    PointSet& operator|=(const Point& a_iv);

    ///
    PointSet& operator|=(const Box& a_box);

    ///and
    PointSet& operator&=(const PointSet& s_sivs);

    ///and
    PointSet& operator&=(const Box& a_box);


    ///not
    PointSet& operator-=(const PointSet& a_sivs);

    ///not
    PointSet& operator-=(const Point& a_iv);

    ///not
    PointSet& operator-=(const Box& a_box);

    ///returns true if
    bool operator==(const PointSet& lhs) const;

    ///
    bool contains(const Point& iv) const;

    ///
    bool contains(const Box& box) const;

    ///
    void grow(int igrow);

    ///
    void grow(int idir, int igrow);

    ///
    void refine(int iref = 2);

    ///
    void coarsen(int iref = 2);

    ///
    void shift(const Point& iv);

    ///
    void clear();

    ///
    Box minBox() const;

    ///
    bool isEmpty() const;

    ///
    void makeEmpty() ;

    ///
    int numPts() const;

  private:

    std::set<Point> m_stdSet;
  };

  ///
  class PointSetIterator
  {
  public:
    ///
    PointSetIterator();

    ///
    PointSetIterator(const PointSet& ivs);

    ///
    void define(const PointSet& ivs);

    ///
    const Point& operator()() const ;

    ///
    bool ok() const;

    ///
    void operator++();

    ///
    void begin();

    ///
    void end();

    ///
    void clear();
  private:
    const PointSet* m_ivs;
    std::set<Point>::iterator m_iter;
  };
}
#include "Proto_PointSetImplem.H"
#endif 
