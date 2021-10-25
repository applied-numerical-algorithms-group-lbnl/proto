using std::shared_ptr;
using std::array;

DisjointBoxLayout::DisjointBoxLayout()
{
}

DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_boxSize);
}

DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const vector<Point>   & a_coarsenedPatches,
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_coarsenedPatches, a_boxSize);
}

void 
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain, 
                          const Point           & a_boxSize)
{
    PROTO_ASSERT(a_problemDomain.box().coarsenable(a_boxSize),
    "DisjointBoxLayout::define | Error: incompatible inputs");    
    Box bxCoarse = a_problemDomain.box().coarsen(a_boxSize);
    BoxIterator bxit(bxCoarse);
    vector<Point> allPoints;
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
    PROTO_ASSERT(a_problemDomain.box().coarsenable(a_boxSize),
    "DisjointBoxLayout::define | Error: incompatible inputs");    
    m_internals = std::make_shared<DBLInternals>();
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
    cout << "end = " << m_end << endl;
}

unsigned int
DisjointBoxLayout::offset(int a_proc) const
{
    return m_internals->startProc[a_proc];
}

///  
// number of boxes in grid (over all procs)
unsigned int
DisjointBoxLayout::size() const
{
    return m_internals->allBoxes.size();
}  

unsigned int
DisjointBoxLayout::localSize() const
{
    return m_internals->startProc[Proto::procID()+1] - m_internals->startProc[Proto::procID()];
}

unsigned int 
DisjointBoxLayout::mySize() const
{
    return localSize();
}

Box
DisjointBoxLayout::operator[](const DataIndex& a_index) const
{
    return box(a_index);
    //return m_internals->allBoxes[a_index.intIndex()].first;
}

Point
DisjointBoxLayout::point(const DataIndex& a_index) const
{
    return m_internals->allBoxes[a_index.intIndex()].first;
}

Box
DisjointBoxLayout::box(const DataIndex& a_index) const
{
    Point p = point(a_index);
    //return Box(p,p).refine(m_boxSize);
    return Box(p*m_boxSize,(p + Point::Ones())*m_boxSize - Point::Ones());
}

int
DisjointBoxLayout::procID(const DataIndex& a_index) const
{
    return m_internals->allBoxes[a_index.intIndex()].second;
}

int
DisjointBoxLayout::procid(const DataIndex& a_index) const
{
    return procID(a_index);
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
DisjointBoxLayout::localIndex(unsigned int a_myIndexInt) const
{
    return DataIndex(m_internals,m_internals->startProc[Proto::procID()] + a_myIndexInt);
}

DataIndex
DisjointBoxLayout::myIndex(unsigned int a_myIndexInt) const
{
    return localIndex(a_myIndexInt);
}

int
DisjointBoxLayout::localIntIndex(const DataIndex& a_di) const
{
    int retval = a_di.intIndex() - m_internals->startProc[Proto::procID()];
    PR_assert(retval >= 0);
    PR_assert(a_di.intIndex() <  m_internals->startProc[Proto::procID()+1]);
    return retval;
}

int
DisjointBoxLayout::myIntIndex(const DataIndex& a_di) const
{
    return localIntIndex(a_di);
}

DataIndex
DisjointBoxLayout::find(const Point& a_pt) const
{
  Point shiftedPt = m_problemDomain.shifted(a_pt);
  // cout << "in find:" << a_pt << " , " << shiftedPt << endl;
  auto val = (m_internals->boxMap).find(mortonIndex(shiftedPt,m_internals->morton));
  if (val==(m_internals->boxMap).end())
        {
          cout << a_pt << " , " << shiftedPt << " , " << m_problemDomain.box() << endl;
        }
    if (val != (m_internals->boxMap).end())
    {
      return DataIndex(m_internals,val->second);
    }
    else 
    {
      return m_end;
    }
}

int
DisjointBoxLayout::maxLocalSize() const
{
    int maxValue = 0;
    for (int ii = 0; ii < Proto::numProc(); ii++)
    {
        int isize = m_internals->startProc[ii+1] - m_internals->startProc[ii];
        maxValue = max(isize, maxValue);
    }
    return maxValue;
}

int
DisjointBoxLayout::minLocalSize() const
{
    int minValue = maxLocalSize();
    for (int ii = 0; ii < Proto::numProc(); ii++)
    {
        int isize = m_internals->startProc[ii+1] - m_internals->startProc[ii];
        minValue = min(isize, minValue);
    }
    return minValue;
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

bool
DisjointBoxLayout::radialSymmetry(Point a_origin) const
{
    auto& boxes = m_internals->allBoxes;
    std::vector<Point> assymmetries; 
    for (auto iter = boxes.begin(); iter != boxes.end(); ++iter)
    {
        Point shifted = (*iter).first - a_origin;
        Point mirrorImage = -(shifted + Point::Ones());
        Point searchPoint = mirrorImage + a_origin;
        if (find(searchPoint) == end())
        {
            assymmetries.push_back((*iter).first);
        }
    }
    if (assymmetries.size() == 0){return true;}
    else {return false;}
}

bool
DisjointBoxLayout::mirrorSymmetry(Point a_origin, int a_coord) const
{
    if (a_coord < 0 || a_coord >= DIM) { return true; }
    
    auto& boxes = m_internals->allBoxes;
    std::vector<Point> assymmetries; 
    for (auto iter = boxes.begin(); iter != boxes.end(); ++iter)
    {
        Point shifted = (*iter).first - a_origin;
        Point mirrorImage = shifted;
        mirrorImage[a_coord] = -(mirrorImage[a_coord] + 1);
        Point searchPoint = mirrorImage + a_origin;
        if (find(searchPoint) == end())
        {
            assymmetries.push_back((*iter).first);
        }
    }
    if (assymmetries.size() == 0){return true;}
    else {return false;}
}

void
DisjointBoxLayout::print(std::string a_name) const
{
    if (Proto::procID() == 0)
    {
        std::cout << "Disjoint Box Layout: " << a_name << std::endl;
        std::cout << std::setw(40) << std::setfill('-') << "-" << std::endl;
        std::cout << std::setfill(' ');
        std::cout << std::setw(10) << std::left << "Domain: "   << m_problemDomain.box() << std::endl;
        std::cout << std::setw(10) << std::left << "Box Size: " << m_boxSize << std::endl;
        std::cout << std::setw(40) << std::setfill('-') << "-" << std::endl;
        std::cout << std::setfill(' ');
        auto boxes = m_internals->allBoxes;
        if (boxes.size() == 0)
        {
            std::cout << "( EMPTY LAYOUT ) " << std::endl;
            std::cout << setfill('\0');
            return;
        }
        std::cout << std::right;
        for (int ii = 0; ii < m_internals->allBoxes.size(); ii++)
        {
            Point p = boxes[ii].first;
            Box b = Box(p,p).refine(m_boxSize);
            std::cout << std::setw(4) << boxes[ii].second << " | ";
            std::cout << b << " | ";
            std::cout << p << std::endl;
        }
        std::cout << setfill(' ') << std::endl;
    }
    
}




