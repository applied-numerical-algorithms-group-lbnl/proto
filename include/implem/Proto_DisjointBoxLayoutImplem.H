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

DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Box             & a_region,
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_region, a_boxSize);
}

DisjointBoxLayout::DisjointBoxLayout(const DisjointBoxLayout& a_input)
{
    m_internals = a_input.m_internals;
    m_problemDomain = a_input.m_problemDomain;
    m_boxSize = a_input.m_boxSize;
    m_end = a_input.m_end;
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
    for (bxit.begin();bxit.ok();++bxit)
    {
        allPoints.push_back(*bxit);
    }
    define(a_problemDomain,allPoints,a_boxSize);
}

void
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain,
                          const Box             & a_region,
                          const Point           & a_boxSize)
{
    PROTO_ASSERT(a_region.coarsenable(a_boxSize),
    "DisjointBoxLayout::define | Error: incompatible inputs");
    Box patchBox = a_region.coarsen(a_boxSize);
    std::vector<Point> patches;
    for (auto iter = patchBox.begin(); iter.ok(); ++iter)
    {
        patches.push_back(*iter);
    }
    define(a_problemDomain, patches, a_boxSize);
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
    //cout << "end = " << m_end << endl;
}

DisjointBoxLayout&
DisjointBoxLayout::operator=(const DisjointBoxLayout& a_input)
{
    m_internals = a_input.m_internals;
    m_problemDomain = a_input.m_problemDomain;
    m_boxSize = a_input.m_boxSize;
    m_end = a_input.m_end;
    return *this;
}

bool
DisjointBoxLayout::operator==(const DisjointBoxLayout& a_input) const
{
    return ((a_input.m_internals == m_internals)&&
            (a_input.m_problemDomain == m_problemDomain)&&
            (a_input.m_boxSize == m_boxSize));
}

const std::vector<std::pair<Point, int>>
DisjointBoxLayout::boxes() const
{
    return m_internals->allBoxes; 
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

Box
DisjointBoxLayout::operator[](const DataIndex& a_index) const
{
    return box(a_index);
}

Point
DisjointBoxLayout::point(const DataIndex& a_index) const
{
    return m_internals->allBoxes[a_index.global()].first;
}

Box
DisjointBoxLayout::box(const DataIndex& a_index) const
{
    Point p = point(a_index);
    return Box(p*m_boxSize,(p + Point::Ones())*m_boxSize - Point::Ones());
}

int
DisjointBoxLayout::procID(const DataIndex& a_index) const
{
    return m_internals->allBoxes[a_index.global()].second;
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
DisjointBoxLayout::find(const Point& a_pt) const
{
    Point periodicImage = a_pt;
    if (patchDomain().image(periodicImage))
    {
        auto val = (m_internals->boxMap).find(mortonIndex(periodicImage, m_internals->morton));
        if (val != (m_internals->boxMap).end())
        {
            return DataIndex(m_internals,val->second);
        }
        else 
        {
            return m_end;
        }
    }
    return m_end;
    
    /*
    Point shiftedPt = patchDomain().shifted(a_pt);
    
    auto val = (m_internals->boxMap).find(mortonIndex(shiftedPt,m_internals->morton));
    if (val != (m_internals->boxMap).end())
    {
        return DataIndex(m_internals,val->second);
    }
    else 
    {
        return m_end;
    }
    */
}

int
DisjointBoxLayout::maxLocalSize() const
{
    int maxValue = 0;
    for (int ii = 0; ii < Proto::numProc(); ii++)
    {
        int isize = m_internals->startProc[ii+1] - m_internals->startProc[ii];
        maxValue = std::max(isize, maxValue);
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
        minValue = std::min(isize, minValue);
    }
    return minValue;
}

bool 
DisjointBoxLayout::coarsenable(const Point& a_refRatio) const
{
    return (((m_problemDomain.box().coarsenable(a_refRatio))) &&
            ((m_boxSize%a_refRatio) == (Point::Zeros())));
}

bool 
DisjointBoxLayout::compatible(const DisjointBoxLayout& a_layout) const
{
    if (m_internals == a_layout.m_internals) {return true; } //same pointer
    return m_internals->compatible(*a_layout.m_internals);
}

bool 
DisjointBoxLayout::compatible(const DataIndex& a_index) const
{
    if (m_internals == a_index.m_internals) { return true; } // same pointer
    return m_internals->compatible(*a_index.m_internals);
}

DisjointBoxLayout 
DisjointBoxLayout::coarsen(const Point& a_refRatio) const
{
    PROTO_ASSERT(coarsenable(a_refRatio),
        "DisjointBoxLayout::coarsen | Error: Invalid coarsening ratio (not coarsenable)");
    DisjointBoxLayout retval;
    retval.m_internals = m_internals;
    retval.m_problemDomain = m_problemDomain.coarsen(a_refRatio);
    retval.m_boxSize = m_boxSize/a_refRatio;
    return retval;
} 

DisjointBoxLayout 
DisjointBoxLayout::coarsen(const Point& a_refRatio, const Point& a_boxSize) const
{
    // The simple case where there are the same number of boxes 
    if (boxSize() / a_refRatio == a_boxSize)
    {
        return coarsen(a_refRatio); 
    }
    MayDay<void>::Error("DisjointBoxLayout::coarsen | Error: Function not implemented.");
    PROTO_ASSERT(coarsenable(a_refRatio),
        "DisjointBoxLayout::coarsen | Error: Invalid coarsening ratio (not coarsenable)");
    auto patches = m_internals->allBoxes;
    for (auto iter = patches.begin(); iter != patches.end(); ++iter)
    {
        Point p = (*iter).first;
        
    }
}

DisjointBoxLayout 
DisjointBoxLayout::refine(const Point& a_refRatio) const
{
    DisjointBoxLayout retval;
    retval.m_internals = m_internals;
    retval.m_problemDomain = m_problemDomain.refine(a_refRatio);
    retval.m_boxSize = m_boxSize*a_refRatio;
    return retval;    
}

Box
DisjointBoxLayout::boundingBox() const
{
    auto& boxes = m_internals->allBoxes;
    if (boxes.size() == 0) {return Box(); }
    Box boundingBox;
    for (auto iter = boxes.begin(); iter != boxes.end(); ++iter)
    {
        Point p = (*iter).first;
        if (boundingBox.empty())
        {
            boundingBox = Box(p, p);
        } else {
            boundingBox &= p;
        }
    }
    return boundingBox.refine(boxSize());
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




