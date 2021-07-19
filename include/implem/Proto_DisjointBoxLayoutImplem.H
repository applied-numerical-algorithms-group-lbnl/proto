using std::shared_ptr;
using std::array;

DisjointBoxLayout::
DisjointBoxLayout()
{
}
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Point           & a_boxSize)
{
  define(a_problemDomain, a_boxSize);
}
void 
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain, 
                          const Point           & a_boxSize)
{
  //cout << "pdBox = " << a_problemDomain.box() << ", a_boxSize = " << a_boxSize << endl;
  PROTO_ASSERT(a_problemDomain.box().coarsenable(a_boxSize), "invalid dbl combo");    
  Box bxCoarse =
    a_problemDomain.box().coarsen(a_boxSize);
  BoxIterator bxit(bxCoarse);
  vector<Point> allPoints;
  //cout << bxCoarse << endl;
  for (bxit.begin();!bxit.end();++bxit)
    {
      allPoints.push_back(*bxit);
    }
  define(a_problemDomain,allPoints,a_boxSize);
}
//
void
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain,
                          const vector<Point>   & a_coarsenedPatches,
                          const Point           & a_boxSize)
{
  PROTO_ASSERT(a_problemDomain.box().coarsenable(a_boxSize), "invalid dbl combo");
  DBLInternals* ptr= new DBLInternals;
  m_internals = shared_ptr<DBLInternals>(ptr);
  m_internals->coarsenedDomain  = a_problemDomain.coarsen(a_boxSize);  
  m_boxSize = a_boxSize;
  m_problemDomain = a_problemDomain;
  BoxMap& boxmap = m_internals->boxMap;
  vector<int>& startproc = m_internals->startProc;
  vector<pair<Point,int> >& allboxes = m_internals->allBoxes;
  vector<Point> patchVec = a_coarsenedPatches;
  mortonSort(patchVec,m_internals->morton);
  int nsegs = numProc();
  int length = patchVec.size() / nsegs;
  int remainder = patchVec.size() % nsegs;
  int ltot = 0;
  for (int k = 0; k < nsegs; k++)
    {
      startproc.push_back(ltot);
      int lengthSeg = length;   
      if (k < remainder)
        {
          lengthSeg = length + 1;
        }
      for (int l = 0;l < lengthSeg;l++)
        {
          boxmap.insert(pair<uint64_t,int>(mortonIndex(patchVec[ltot],m_internals->morton),ltot));
          auto pr = pair<Point,int>(patchVec[ltot],k);
          allboxes.push_back(pr);
          ltot++;
        }
    }
  // Adding one more element to startproc avoid a conditional in computing the number of boxes 
  // assigned to each processor.
  startproc.push_back(patchVec.size());
  m_end = DataIndex(m_internals,patchVec.size());
}
///  
// number of boxes in grid (over all procs)
unsigned int
DisjointBoxLayout::size() const
{
  return m_internals->allBoxes.size();
}  
unsigned int 
DisjointBoxLayout::mySize() const
{
  return m_internals->startProc[procID()+1] - m_internals->startProc[procID()];
}
Point
DisjointBoxLayout::operator[](const DataIndex& a_index) const
{
  return m_internals->allBoxes[a_index.intIndex()].first;
}
int
DisjointBoxLayout::procid(const DataIndex& a_index) const
{
  return m_internals->allBoxes[a_index.intIndex()].second;
}
DataIndex
DisjointBoxLayout::index(Point& a_coarsenedPt) const
{
  Point shiftedPt = m_internals->coarsenedDomain.shifted(a_coarsenedPt);
  return DataIndex(m_internals,m_internals->boxMap[mortonIndex(shiftedPt,m_internals->morton)]);
}
DataIndex
DisjointBoxLayout::index(unsigned int a_indexInt) const
{
  return DataIndex(m_internals,a_indexInt);
}
DataIndex
DisjointBoxLayout::myIndex(unsigned int a_myIndexInt) const
{
  return DataIndex(m_internals,m_internals->startProc[procID()] + a_myIndexInt);
}
int
DisjointBoxLayout::myIntIndex(const DataIndex& a_di) const
{
  int retval = a_di.intIndex() - m_internals->startProc[procID()];
  //cout << "in myIntIndex:" << a_di.intIndex() << " , " << m_internals->startProc[procID()+1] << endl;
  PR_assert(retval >= 0); PR_assert(a_di.intIndex() <  m_internals->startProc[procID()+1]);
  return retval;
}
DataIndex
DisjointBoxLayout::find(const Point& a_pt) const
{
  auto shiftedPt = m_internals->coarsenedDomain.shifted(a_pt);
  auto val = (m_internals->boxMap).find(mortonIndex(shiftedPt,m_internals->morton));
  if (val != (m_internals->boxMap).end())
    {
      return DataIndex(m_internals,val->second);
    }
  else
    {
      return this->end();
    }
}
bool 
DisjointBoxLayout::coarsenable(const Point& a_refrat) const
{
  return (((m_problemDomain.box().coarsenable(a_refrat))) &&
           ((m_boxSize%a_refrat) == (Point::Zeros())));
}
DisjointBoxLayout 
DisjointBoxLayout::coarsen(const Point& a_refrat) const
{
  PROTO_ASSERT(coarsenable(a_refrat),"DBL not coarsenable");
  DisjointBoxLayout retval;
  retval.m_internals = m_internals;
  retval.m_problemDomain = m_problemDomain.coarsen(a_refrat);
  retval.m_boxSize = m_boxSize/a_refrat;
  return retval;
} 
DisjointBoxLayout 
DisjointBoxLayout::refine(const Point& a_refrat) const
{
  DisjointBoxLayout retval;
  retval.m_internals = m_internals;
  retval.m_problemDomain = m_problemDomain.refine(a_refrat);
  retval.m_boxSize = m_boxSize*a_refrat;
  return retval;    
}
