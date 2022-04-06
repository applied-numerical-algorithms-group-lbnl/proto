using std::shared_ptr;
using std::array;

// Default Constructor
DisjointBoxLayout::DisjointBoxLayout()
{
}

// Simple Constructor
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_boxSize);
}

// General Constructor
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const vector<Point>   & a_coarsenedPatches,
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_coarsenedPatches, a_boxSize);
}

// Sub-Region Constructor
DisjointBoxLayout::DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                                     const Box             & a_region,
                                     const Point           & a_boxSize)
{
    define(a_problemDomain, a_region, a_boxSize);
}

// Copy Constructor
DisjointBoxLayout::DisjointBoxLayout(const DisjointBoxLayout& a_input)
{
    m_internals = a_input.m_internals;
    m_problemDomain = a_input.m_problemDomain;
    m_boxSize = a_input.m_boxSize;
    m_end = a_input.m_end;
}

// Define (Simple)
void 
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain, 
                          const Point           & a_boxSize)
{
    Box bxCoarse = a_problemDomain.box().coarsen(a_boxSize);
    BoxIterator bxit(bxCoarse);
    vector<Point> allPoints;
    for (bxit.begin();bxit.ok();++bxit)
    {
        allPoints.push_back(*bxit);
    }
    define(a_problemDomain,allPoints,a_boxSize);
}

// Define (Sub-Domain)
void
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain,
                          const Box             & a_region,
                          const Point           & a_boxSize)
{
    PROTO_ASSERT(a_problemDomain.box().contains(a_region),
        "DisjointBoxLayout::define | Error: Sub-region is not a subset of the domain.");
    Box patchBox = a_region.coarsen(a_boxSize);
    std::vector<Point> patches;
    for (auto iter = patchBox.begin(); iter.ok(); ++iter)
    {
        patches.push_back(*iter);
    }
    define(a_problemDomain, patches, a_boxSize);
}

// Define (General)
void
DisjointBoxLayout::define(const ProblemDomain   & a_problemDomain,
                          const vector<Point>   & a_coarsenedPatches,
                          const Point           & a_boxSize)
{
    PROTO_ASSERT(a_problemDomain.box().coarsenable(a_boxSize),
        "DisjointBoxLayout::define | Error: Sub-region cannot be tiled by input box size.");    
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
}

// Assignment Opoerator
DisjointBoxLayout&
DisjointBoxLayout::operator=(const DisjointBoxLayout& a_input)
{
    m_internals = a_input.m_internals;
    m_problemDomain = a_input.m_problemDomain;
    m_boxSize = a_input.m_boxSize;
    m_end = a_input.m_end;
    return *this;
}

// Equality Operator
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
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::operator[] | Error: \
        Input DataIndex is not compatible with this layout.");
    return box(a_index);
}

Point
DisjointBoxLayout::point(const DataIndex& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::point | Error: Input DataIndex is not compatible with this layout.");
    return m_internals->allBoxes[a_index.global()].first;
}

Box
DisjointBoxLayout::box(const DataIndex& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::box | Error: Input DataIndex is not compatible with this layout.");
    Point p = point(a_index);
    return Box(p*m_boxSize,(p + Point::Ones())*m_boxSize - Point::Ones());
}

int
DisjointBoxLayout::procID(const DataIndex& a_index) const
{
    PROTO_ASSERT(compatible(a_index),
        "DisjointBoxLayout::procID | Error: Input DataIndex is not compatible with this layout.");
    return m_internals->allBoxes[a_index.global()].second;
}

DataIndex
DisjointBoxLayout::index(Point& a_patchPoint) const
{
    PROTO_ASSERT(patchDomain().contains(a_patchPoint),
        "DisjointBoxLayout::index | Error: \
        Point does not represent a valid Patch in the layout or any periodic image.");
    Point shiftedPt = m_internals->coarsenedDomain.image(a_patchPoint);
    return DataIndex(m_internals,m_internals->boxMap[mortonIndex(shiftedPt,m_internals->morton)]);
}

DataIndex
DisjointBoxLayout::index(unsigned int a_indexInt) const
{
    PROTO_ASSERT(a_indexInt < size(),
        "DisjointBoxLayout::index | Error: \
        Global integer index %u is out of bounds.", a_indexInt);
    return DataIndex(m_internals,a_indexInt);
}

DataIndex
DisjointBoxLayout::localIndex(unsigned int a_indexInt) const
{
    PROTO_ASSERT(a_indexInt < localSize(),
        "DisjointBoxLayout::index | Error: \
        Global integer index %u is out of bounds.", a_indexInt);
    return DataIndex(m_internals,m_internals->startProc[Proto::procID()] + a_indexInt);
}

DataIndex
DisjointBoxLayout::find(const Point& a_pt) const
{
    if (!patchDomain().contains(a_pt)) { return m_end; }
    else
    {
        Point periodicImage = patchDomain().image(a_pt);
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
}

bool
DisjointBoxLayout::contains(const Point& a_pt) const
{
    return find(a_pt) != end();
}

bool
DisjointBoxLayout::onDomainBoundary(const Point& a_pt) const
{
    PROTO_ASSERT(find(a_pt) != end(),
        "DisjointBoxLayout::onDomainBoundary | Error: \
        Point is not in this layout.");
    Box K = Box::Kernel(1).shift(a_pt);
    
    for (auto biter = K.begin(); biter != K.end(); ++biter)
    {
        if (*biter == a_pt) {continue;}
        if (!patchDomain().contains(*biter)) { return true; }
    }
    return false;
}

bool
DisjointBoxLayout::onLevelBoundary(const Point& a_pt) const
{
    PROTO_ASSERT(find(a_pt) != end(),
        "DisjointBoxLayout::onDomainBoundary | Error: \
        Point is not in this layout.");
    Box K = Box::Kernel(1).shift(a_pt);
    
    for (auto biter = K.begin(); biter != K.end(); ++biter)
    {
        if (*biter == a_pt) {continue;}
        if (find(*biter) == end()) { return true; }
    }
    return false;
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
    return DisjointBoxLayout();
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
    using Proto::pout;

    pout() << "Disjoint Box Layout: " << a_name << std::endl;
    pout() << std::setw(40) << std::setfill('-') << "-" << std::endl;
    pout() << std::setfill(' ');
    pout() << std::setw(10) << std::left << "Domain: "   << m_problemDomain.box() << std::endl;
    pout() << std::setw(10) << std::left << "Box Size: " << m_boxSize << std::endl;
    pout() << std::setw(40) << std::setfill('-') << "-" << std::endl;
    pout() << std::setfill(' ');
    auto boxes = m_internals->allBoxes;
    if (boxes.size() == 0)
    {
        pout() << "( EMPTY LAYOUT ) " << std::endl;
        pout() << setfill('\0');
        return;
    }
    pout() << std::right;
    for (int ii = 0; ii < m_internals->allBoxes.size(); ii++)
    {
        Point p = boxes[ii].first;
        Box b = Box(p,p).refine(m_boxSize);
        pout() << std::setw(4) << boxes[ii].second << " | ";
        pout() << b << " | ";
        pout() << p << std::endl;
    }
    pout() << setfill(' ') << std::endl;
}

inline std::ostream& operator<< (std::ostream& os, const DisjointBoxLayout& a_dbl)
{
    os << "DisjointBoxLayout: \n" ;
    os << a_dbl.domain() << "\n";
    auto internals = a_dbl.m_internals;
    os << "coarsened Domain = " << internals->coarsenedDomain << "\n" ;
    os << "Points and ranks: \n";
    for (int k = 0; k < internals->allBoxes.size();k++)
    {
        os << internals->allBoxes[k].first 
            << " , " << internals->allBoxes[k].second << "\n";
    }
    os << "StartProc vector: \n";
    for (int k = 0; k < internals->startProc.size() ; k++)
    {
        os << "Starting patch number for processor " << k 
            << " = " << internals->startProc[k] << "\n";
    }
    return os;
}



