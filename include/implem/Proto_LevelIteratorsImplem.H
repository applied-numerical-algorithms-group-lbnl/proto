//
namespace Proto
{
    NeighborIterator::NeighborIterator(const DisjointBoxLayout& a_dbl,
                                     const Box& a_bx)
  {
    m_dbl = a_dbl;
    m_bx = a_bx;
    m_current = DataIndex(a_dbl,0);
    m_end.define(a_dbl);
    
    // Set up BoxIterator.
    ProblemDomain dProblem = m_dbl.problemDomain();
    Point boxSize = m_dbl.boxSize();
    ProblemDomain dProbCoarse = dProblem.coarsen(boxSize);
    Box bxPoints = dProbCoarse&m_bx.coarsen(m_dbl.boxSize());
    m_bxit = BoxIterator(bxPoints);
  }
  NeighborIterator&
  NeighborIterator::begin()
  {
    m_bxit.begin();
    while ((m_dbl.find(bxit()) == m_dbl.end())&&(!m_bxit.end()))
      {
        ++m_bxit;
      }
    if (!m_bxit.end())
      {
        Point shiftCoarse = dProbCoarse.shift(bxit());
        m_shift = shiftCoarse*m_dbl.boxSize();
        m_current = m_dbl.find(m_bxit());
      }
    else
      {
        m_current = m_end;
      }
    return *this;
  }
  NeighborIterator&
  NeighborIterator::operator++()
  {
    ++m_bxit;
    while ((m_dbl.find(m_bxit()) == m_dbl.end())&&(!m_bxit.end()))
      {
        ++m_bxit;
      }
    if (!m_bxit.end())
      {
        Point shiftCoarse = dProbCoarse.shift(bxit());
        m_shift = shiftCoarse*m_dbl.boxSize();
        m_current = m_dbl.find(m_bxit());
      }
    else
      {
        m_current = m_end;
      }
    return *this;
  }
  Box
  NeighborIterator::box();
  {
    Box bx(m_bxit(),m_bxit());
    return bx.refine(m_dbl.BoxSize);
  }
  Point
  NeighborIterator::neighborShift()
  {
    return m_shift;
  }
  DataIterator::DataIterator(DisjointBoxLayout a_dbl)
  {
    m_dblInternals = a_dbl.m_internals;
    m_currentIndex = DataIndex(a_dbl);
    m_end.defineEnd(a_dbl);
  }
  DataIterator&
  DataIterator::begin()
  {
    m_currentIndex = m_dbl.myIndex(0);
    return *this;
  }
  DataIterator&
  DataIterator::operator++
  {
    currentInt = m_current.currentInt();
    if (currentInt < m_dbl.mySize())
      {
        currentIndex = m_dbl.index(m_currentInt);
      }
    else
      {
        m_currentIndex = m_end;
      }
    return *this;
  }
  DataIndex
  DataIterator::operator*
  {
    return m_currentIndex;
  } 
  Point
  DataIterator::operator()
  {
    return m_dbl[m_currentIndex];
  }
  Box
  DataIterator::box()
  {
    Point lowCorner = m_dbl[m_currentIndex]*m_dbl.boxSize();
    return Box(lowCorner,lowCorner+m_dbl.boxSize() - Point::Ones());
  }
#if 0
  FineFRIterator::FineFRIterator(
                                 const DataIndex& a_fineDit,
                                 const Point& a_refRatio)
  {
  }
  FRFineIterator&
  FRFineIterator::begin()
  {
  }
  DataIndex&
  FRFineIterator::end()
  {
  }
  FRFineIterator&
  FRFineIterator::operator++()
  {
  }
  DataIndex&
  FRFineIterator::operator()
  {
  }
  Box
  FRFineIterator::box()
  {
  }
  Point
  FRFineIterator::srcShift()
  {
  }
  CoarseFRIterator::CoarseFRIterator(const DataIndex& a_coarseDit,
                                     const Point& a_refRatio)
  {
  }
  CoarseFRIterator&
  CoarseFRIterator::begin()
  {
  }
  CoarseFRIterator&
  CoarseFRIterator::end()
  {
  }
  CoarseFRIterator&
  CoarseFRIterator::operator++()
  {
  }
  DataIndex
  CoarseFRIterator::operator()
  {
  }
  Box
  CoarseFRIterator::box()
  {
  }
  Point
  CoarseFRIterator::srcShift()
  {
  }
  DataIndex
  CoarseFRIterator::fromIndex()
  {
  }
#endif
}
