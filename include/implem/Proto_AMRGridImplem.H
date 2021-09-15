
void
AMRGrid::buffer(LevelTagData& a_tags, int a_bufferSize)
{
    AMRGrid::buffer(a_tags, Point::Ones(a_bufferSize));
}

void
AMRGrid::buffer(LevelTagData& a_tags, Point a_bufferSize)
{
    Box K(-a_bufferSize, a_bufferSize);
    Stencil<char> S;
    char w = 1;
    for (auto iter = K.begin(); !iter.end(); ++iter)
    {
        S += w*Shift(*iter);
    }
    a_tags.exchange();
    for (auto iter = a_tags.begin(); iter.ok(); ++iter)
    {
        auto& patch = a_tags[*iter];
        TagData tmp(patch.box(), false); //not using the Proto::Stack for this.
        patch.copyTo(tmp);
        patch |= S(tmp);
        
        forallInPlace(
        [=] PROTO_LAMBDA (Var<char>& a_data)
        {
            a_data(0) = min((int)a_data(0), 1);
        }, patch);
    }
}

void
AMRGrid::define(std::vector<DisjointBoxLayout>& a_layouts, int a_maxLevels)
{
    PROTO_ASSERT(a_maxLevels > 0,
        "AMRGrid::define | Error: maxLevels should not be zero.");
    PROTO_ASSERT(a_layouts.size() <= a_maxLevels,
        "AMRGrid::define | Error: Input layout array has %lu levels, but maxLevels is set to %i.",
        a_layouts.size(), a_maxLevels);
    m_maxLevels = a_maxLevels;
    m_layouts.clear();
    for (int ii = 0; ii < a_layouts.size(); ii++)
    {
        m_layouts.push_back(a_layouts[ii]); //copying here is fine; dbl internals do their job
    }
    m_defined = true;
}

void
AMRGrid::define(DisjointBoxLayout& a_layout, int a_maxLevels)
{
    PROTO_ASSERT(a_maxLevels > 0,
        "AMRGrid::define | Error: maxLevels should not be zero.");
    m_maxLevels = a_maxLevels;
    std::vector<DisjointBoxLayout> layouts;
    layouts.push_back(a_layout);
    define(layouts, a_maxLevels);
}

const DisjointBoxLayout&
AMRGrid::operator[](unsigned int a_level) const
{
    PROTO_ASSERT((a_level >= 0) && (a_level < m_layouts.size()),
        "AMRGrid::operator[] | Error: level index %i out of bounds. Should be in [0, %lu)",
        a_level, m_layouts.size());

    return m_layouts[a_level];
}

DisjointBoxLayout&
AMRGrid::operator[](unsigned int a_level)
{
    PROTO_ASSERT((a_level >= 0) && (a_level < m_layouts.size()),
        "AMRGrid::operator[] | Error: level index %i out of bounds. Should be in [0, %lu)",
        a_level, m_layouts.size());

    return m_layouts[a_level];
}

void
AMRGrid::regrid(LevelTagData& a_bufferedTags, unsigned int a_level)
{
    PROTO_ASSERT(a_level < m_layouts.size(),
        "AMRGrid::regrid | Error: Input level %u exceeds size %lu of AMRGrid.", a_level, m_layouts.size());
    //TODO: Maybe regridding on the finest level should be a null-op instead of an error.  
    PROTO_ASSERT(a_level+1 < m_maxLevels,
        "AMRGrid::regrid | Error: Attempt to regrid on the finest level.");
    Point boxSize = operator[](a_level).boxSize();
    regrid(a_bufferedTags, a_level, boxSize);

}
void
AMRGrid::regrid(LevelTagData& a_bufferedTags, unsigned int a_level, Point a_fineBoxSize)
{
    PROTO_ASSERT(a_level < m_layouts.size(),
        "AMRGrid::regrid | Error: Input level %u exceeds size %lu of AMRGrid.", a_level, m_layouts.size());
    //TODO: Maybe regridding on the finest level should be a null-op instead of an error.  
    PROTO_ASSERT(a_level+1 < m_maxLevels,
        "AMRGrid::regrid | Error: Attempt to regrid on the finest level.");
    auto& crseLayout = operator[](a_level);
    PROTO_ASSERT(crseLayout == a_bufferedTags.layout(),
        "AMRGrid::regrid | Error: Tag data layout does not match AMRGrid layout on level %u.", a_level);
   
    Point patchRefRatio = crseLayout.boxSize() / a_fineBoxSize * PR_AMR_REFRATIO;
    Box crseBitDomain = crseLayout.pointDomain();
    Box fineBitDomain = crseBitDomain.refine(patchRefRatio);

    FinitePointSet taggedPatches(fineBitDomain, crseLayout.domain().periodicflags());
    for ( auto iter = a_bufferedTags.begin(); iter.ok(); ++iter)
    {
        auto& patch = a_bufferedTags[*iter];
        Point tagRefRatio = crseLayout.boxSize() / patchRefRatio;
        TagData tags(iter.box().coarsen(tagRefRatio));
        tags |= Stencil<char>::AvgDown(tagRefRatio)(patch);
        
        for (auto biter = tags.box().begin(); biter.ok(); ++biter)
        {
            if (tags(*biter) != 0)
            {
                taggedPatches.add(*biter);
            }
        }
    }
  
    ProblemDomain fineProblemDomain = crseLayout.domain().refine(Point::Ones(PR_AMR_REFRATIO));
    DisjointBoxLayout fineLayout(fineProblemDomain, taggedPatches.points(), a_fineBoxSize);
    if (a_level + 1 == m_layouts.size())
    {
        // adding a new level of refinement with the next coarser level's boxSize
        m_layouts.push_back(fineLayout);
    } else {
        // regridding an existing level of refinement
        m_layouts[a_level + 1] = fineLayout;
    }
}

void
AMRGrid::enforceNesting(unsigned int a_level, int a_nestingDistance)
{
    // input level is the "middle" level (e.g. the one being modified)
    PROTO_ASSERT(a_level+1 < m_layouts.size(),
        "AMRGrid::enforceNesting | Error: Invalid input level %u", a_level);
    PROTO_ASSERT(a_level > 0,
        "AMRGrid::enforceNesting | Error: Invalid input level %u", a_level);

    auto& layout = m_layouts[a_level];
    auto& crseLayout = m_layouts[a_level-1];
    auto& fineLayout = m_layouts[a_level+1];
    
    ProblemDomain patchDomain = layout.domain().coarsen(layout.boxSize());

    const auto& patches = layout.dblInternals()->allBoxes;
    const auto& finePatches = fineLayout.dblInternals()->allBoxes;
    
    std::set<Point> points; 
    
    // patch-space refinement ratio between this level and fine level
    Point finePatchRatio = layout.boxSize() / fineLayout.boxSize() * PR_AMR_REFRATIO;

    // Add this layout's patches in the coarse point space
    for (auto iter = patches.begin(); iter != patches.end(); ++iter)
    {
        points.insert((*iter).first);
    }
   
    // Union with the fine layout's patches in the coarse point space 
    int nestingVect[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        int n = ceil((double)(a_nestingDistance * PR_AMR_REFRATIO) / (double)(fineLayout.boxSize()[ii]));
        nestingVect[ii] = n;
    }
    Point growRadius(nestingVect);
    for (auto iter = finePatches.begin(); iter != finePatches.end(); ++iter)
    {
        Point finePoint = (*iter).first;
        auto K = Box(-growRadius, growRadius);
        K = K.shift(finePoint);
        for (auto biter = K.begin(); biter.ok(); ++biter)
        {
            Point pi = (*biter) / finePatchRatio;
            if (patchDomain.image(pi))
            {
                points.insert(pi);
            }
        }
    }

    std::vector<Point> pointVect;
    for (auto iter = points.begin(); iter != points.end(); ++iter)
    {
        pointVect.push_back(*iter);
    }
    
    // Redefine the layout 
    layout.define(layout.domain(), pointVect, layout.boxSize());
}


