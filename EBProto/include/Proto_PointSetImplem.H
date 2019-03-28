#ifndef __Proto_PointSetImplem__
#define __Proto_PointSetImplem__
namespace Proto
{
///
PointSet::
PointSet(const Box& a_box)
{
  for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
  {
    *this |= *bit;
  }
}
///
  PointSet::
  PointSet(const PointSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
  }
///
  void 
  PointSet::
  define(const Box& a_box)
  {
    *this = PointSet(a_box);
  }
///
  void 
  PointSet::
  define(const PointSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
  }
///
  PointSet& 
  PointSet::
  operator=(const PointSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
    return *this;
  }
///
  PointSet& 
  PointSet::
  operator|=(const PointSet& a_sivs)
  {
    const std::set<IntVect, lex_compare_iv> inputset = a_sivs.m_stdSet;
    std::set<IntVect,lex_compare_iv>::iterator it;
    for(it = inputset.begin(); it!=  inputset.end(); ++it)
    {
      m_stdSet.insert(*it);
    }
    return *this;
  }
///
  PointSet & 
  PointSet::
  operator|=(const IntVect& a_iv)
  {
    m_stdSet.insert(a_iv);
    return *this;
  }
///
  PointSet& 
  PointSet::
  operator|=(const Box& a_box)
  {
    for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
    {
      m_stdSet.insert(*bit);
    }
    return *this;
  }
///
  PointSet& 
  PointSet::
  operator&=(const PointSet& a_sivs)
  {
    if(&a_sivs != this)
    {
      std::set<IntVect, lex_compare_iv> newSet;
      std::set<IntVect, lex_compare_iv>::iterator it;
      for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        if(contains(iv) && a_sivs.contains(iv))
        {
          newSet.insert(iv);
        }
      }
      m_stdSet = newSet;
    }
    return *this;
  }

///and
  PointSet& 
  PointSet::
  operator&=(const Box& a_box)
  {
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      if(!a_box.contains(iv))
      {
        m_stdSet.erase(it);
      }
    }
    return *this;
  }
///and
  PointSet& 
  PointSet::
  operator&=(const ProblemDomain& a_domain)
  {
    return (*this &= a_domain.domainBox());
  }
///not
  PointSet& 
  PointSet::
  operator-=(const PointSet& a_sivs)
  {
    std::set<IntVect, lex_compare_iv>::iterator it;
    //leaving out the ++it because  erase 
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); )
    {
      if(a_sivs.contains(*it))
      {
        m_stdSet.erase(it++);
      }
      else
      {
        ++it;
      }
    }
    return *this;
  }
///not
  PointSet& 
  PointSet::
  operator-=(const IntVect& a_iv)
  {
    if(contains(a_iv))
    {
      m_stdSet.erase(m_stdSet.find(a_iv));
    }
    return *this;
  }
///not
  PointSet& 
  PointSet::
  operator-=(const Box& a_box)
  {
    for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
    {
      *this -= *bit;
    }
    return *this;
  }
///
  bool 
  PointSet::
  operator==(const PointSet& a_lhs) const
  {
    if(a_lhs.m_stdSet.size() != m_stdSet.size())
    {
      return false;
    }

    bool retval = true;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
    {
      if((!contains(*it)) || (!a_lhs.contains(*it)))
      {
        retval = false;
        break;
      }
    }
    return retval;
  }

///
  bool 
  PointSet::
  contains(const IntVect& a_iv) const
  {
    std::set<IntVect, lex_compare_iv>::iterator it = m_stdSet.find(a_iv);
    return (it != m_stdSet.end());
  }

///
  bool 
  PointSet::
  contains(const Box& a_box) const
  {
    bool retval = true;
    for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
    {
      if(!contains(*bit))
      {
        retval = false;
        break;
      }
    }
    return retval;
  }

///
  void 
  PointSet::
  grow(int igrow)
  {
    PointSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.grow(igrow);
      newSet |= grid;
    }
    *this = newSet;
  }

///
  void 
  PointSet::
  grow(int idir, int igrow)
  {
    PointSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.grow(idir, igrow);
      newSet |= grid;
    }
    *this = newSet;
  }

///
  void 
  PointSet::
  growHi()
  {
    for(int idir = 0; idir < SpaceDim; idir++)
    {
      growHi(idir);
    }
  }

///
  void 
  PointSet::
  growHi(int a_dir)
  {
    PointSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.growHi(a_dir);
      newSet |= grid;
    }
    *this = newSet;
  }

///
  void 
  PointSet::
  refine(int iref)
  {
    PointSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.refine(iref);
      newSet |= grid;
    }
    *this = newSet;
  }

///
  void 
  PointSet::
  coarsen(int iref)
  {
    std::set<IntVect, lex_compare_iv> newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      Box grid(iv, iv);
      grid.coarsen(iref);
      for(BoxIterator bit = a_box.begin(); bit != a_box.end();  ++bit)
      {
        newSet.insert(*bit);
      }
    }
    m_stdSet = newSet;
  }

///
  void 
  PointSet::
  shift(const IntVect& a_iv)
  {
    std::set<IntVect, lex_compare_iv> newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      IntVect iv = *it;
      iv.shift(a_iv);
      newSet.insert(iv);
    }
    m_stdSet = newSet;
  }

///
  void 
  PointSet::
  clear()
  {
    std::set<IntVect, lex_compare_iv> newSet;
    m_stdSet = newSet;
  }

///
  Box 
  PointSet::
  minBox() const
  {
    int bignum = 100000;
    IntVect lo = bignum*IntVect::Unit;
    IntVect hi =-bignum*IntVect::Unit;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      const IntVect& iv = *it;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        lo[idir] = Min(lo[idir], iv[idir]);
        hi[idir] = Max(hi[idir], iv[idir]);
      }
    }
  
    Box retval(lo, hi);
    return retval;
  }

///
  bool 
  PointSet::
  isEmpty() const
  {
    return (m_stdSet.size() == 0);
  }
///
  void 
  PointSet::
  getVectorIV(Vector<IntVect>& a_vect) const
  {
    a_vect.resize(m_stdSet.size());
    int bignum = 100000;
    IntVect lo = bignum*IntVect::Unit;
    IntVect hi =-bignum*IntVect::Unit;
    std::set<IntVect, lex_compare_iv>::iterator it;
    int ivec = 0;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
    {
      a_vect[ivec] = *it;
      ivec++;
    }
  }
///
  void 
  PointSet::
  makeEmpty() 
  {
    clear();
  }

///
  int
  PointSet::
  numPts() const
  {
    return m_stdSet.size();
  }
///
  int 
  PointSet::
  linearSize() const
  {
    Vector<IntVect> vect;
    getVectorIV(vect);
    int retval = linearListSize(vect);
    return retval;
  }

///
  void 
  PointSet::
  define(const Vector<IntVect>& a_vect)
  {
    makeEmpty();
    for(int ivec = 0; ivec  < a_vect.size(); ivec++)
    {
      (*this ) |= a_vect[ivec];
    }
  }

///
  void 
  PointSet::
  linearIn(const void* const inBuf)
  {
    Vector<IntVect> vect;
    linearListIn(vect, inBuf);
    define(vect);
  }

///
  void 
  PointSet::
  linearOut(void* const a_outBuf) const
  {
    Vector<IntVect> vect;
    getVectorIV(vect);
    linearListOut(a_outBuf, vect);
  }


  PointSetIterator::
  PointSetIterator()
  {
    m_ivs = NULL;
  }

///
  PointSetIterator::
  PointSetIterator(const PointSet& ivs)
  {
    m_ivs = &ivs;
    m_iter = m_ivs->m_stdSet.begin();
  }

///
  void 
  PointSetIterator::
  define(const PointSet& a_ivs)
  {
    m_ivs = &a_ivs;
    m_iter = m_ivs->m_stdSet.begin();
  }

///
  const IntVect& 
  PointSetIterator::
  operator()() const 
  {
    return *m_iter;
  }

///
  bool 
  PointSetIterator::
  ok() const
  {
    PROTO_ASSERT(m_ivs != NULL);
    return (m_iter != m_ivs->m_stdSet.end());
  }

///
  void 
  PointSetIterator::
  operator++()
  {
    m_iter++;
  }

///
  void 
  PointSetIterator::
  begin()
  {
    PROTO_ASSERT(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.begin();
  }

///
  void 
  PointSetIterator::
  end()
  {
    PROTO_ASSERT(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.end();
  }

///
  void 
  PointSetIterator::
  clear()
  {
    m_ivs = NULL;
  }
}
#endif
