
void
AMRGrid::buffer(LevelTagData& a_tags, int a_bufferSize)
{
    Box K = Box::Kernel(a_bufferSize);
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
        TagData tmp(patch.box());
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
AMRGrid::define(std::vector<DisjointBoxLayout>& a_layouts)
{
    m_layouts.clear();
    for (int ii = 0; ii < m_layouts.size(); ii++)
    {
        m_layouts.push_back(a_layouts[ii]); //copying here is fine; dbl internals do their job
    }
    m_defined = true;
}

void
AMRGrid::define(DisjointBoxLayout& a_layout)
{
    std::vector<DisjointBoxLayout> layouts;
    layouts.push_back(a_layout);
    define(layouts);
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
    auto& crseLayout = operator[](a_level);
    PROTO_ASSERT(crseLayout == a_bufferedTags.layout(),
        "AMRGrid::regrid | Error: Tag data layout does not match AMRGrid layout on level %u.", a_level);
    
    
    //TODO: Make sure a_level is one less than the max level or lower.
    
    std::vector<Point> taggedPatches;
    for ( auto iter = a_bufferedTags.begin(); iter.ok(); ++iter)
    {
        auto& patch = a_bufferedTags[*iter];
        auto tagged = patch.absMax();
        if (tagged > 0)
        {
            Box b = crseLayout[*iter];
            Point crseBoxSize = crseLayout.boxSize();
            b = b.coarsen(crseBoxSize).refine(PR_AMR_REFRATIO);
        
            for (auto biter = b.begin(); biter.ok(); ++biter)
            {
                taggedPatches.push_back(*biter);
            }
        }
    }
  
    ProblemDomain problemDomain(crseLayout.domain().box().refine(PR_AMR_REFRATIO), true);
    DisjointBoxLayout fineLayout(problemDomain, taggedPatches, crseLayout.boxSize());
    
    if (a_level + 1 == m_layouts.size())
    {
        // adding a new level of refinement
        m_layouts.push_back(fineLayout);
    } else {
        // regridding an existing level of refinement
        m_layouts[a_level + 1] = fineLayout;
    }

    //TODO: Enforce proper nesting of lower level grids. 
}








