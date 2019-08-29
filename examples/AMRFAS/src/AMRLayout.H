#pragma once

#define AMR_REFRATIO 2
#define MG_REFRATIO 2

#define PRE_RELAX 4
#define POST_RELAX 4
#define BOTTOM_RELAX 20

#define MAXBOXSIZE 32
#define BLOCKFACTOR 2

#include "Proto.H"
#include "ProtoInterface.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "AMRIO.H"

/// Build Layout
/**
    Shortcut utility function to build a Chombo::DisjointBoxLayout with input
    domain Box, periodicity, and optional load balance / domain splitting.

    \param a_layout     An existing DisjointBoxLayout. Usually this is default constructed.
    \param a_box        A Chombo::Box (or Proto::Box) describing the domain of the layout.
    \param a_periodic   A Proto::Point whose ith component is 0 for non-periodic bcs
    \param a_doSplit    Defaults to true in which case the domain is split according to
                        MAXBOXSIZE and BLOCKFACTOR
*/
void buildLayout(DisjointBoxLayout& a_layout,
        const Box& a_box,
        const Proto::Point& a_periodic = Proto::Point::Ones(),
        bool  a_doSplit = true)
{
    PR_TIME("::buildLayout (DBL construction utility found in AMRData.H)");
    bool periodic[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        periodic[ii] = (bool)a_periodic[ii];
    }
    ProblemDomain domain(a_box, periodic);
    Vector<Box> boxVect;
    Vector<int> procAssign;
    if (a_doSplit)
    {
        domainSplit(domain, boxVect, MAXBOXSIZE, BLOCKFACTOR);

        procAssign.resize(boxVect.size(), 0);

        LoadBalance(procAssign, boxVect);
    } else {
        boxVect.push_back(a_box);
        procAssign.push_back(0);
    }
    a_layout.define(boxVect, procAssign, domain);
}

/// Coarsen Layout
/**
    Initialize a DisjointBoxLayout as a coarsened version of an existing layout.
    Always does domain splitting (e.g. coarse and fine layouts are not guaranteed to have
    compatable dataIterators)

    \param a_coarseLayout   An uninitialized DisjointBoxLayout
    \param a_fineLayout     An existing, defined DisjointBoxLayout
    \param a_refRatio       Ratio by which to coarsen
*/
void coarsenGrid(DisjointBoxLayout& a_coarseLayout,
                 const DisjointBoxLayout& a_fineLayout,
                 unsigned int a_refRatio)
{
    PR_TIME("::coarsenGrid (DBL utility found in AMRData.H)");
    ProblemDomain domain = a_fineLayout.physDomain();
    domain.coarsen(a_refRatio);

    Vector<Box> boxVect;
    domainSplit(domain, boxVect, MAXBOXSIZE, BLOCKFACTOR);
    Vector<int> procAssign(boxVect.size(), 0);
    LoadBalance(procAssign, boxVect);
    a_coarseLayout.define(boxVect, procAssign, domain);
}

/// Adaptive Mesh Refinement Layout
/**
    Encapsulates a hierarchy of unions of rectangular domains (e.g. an AMR hierarchy).
    Currently this hierarchy is non-branching (e.g. 1 rectangular domain at each level)
    however proper branching will be implemented in the near future. 
*/

class AMRLayout {
public:
    /// Default Constructor
    AMRLayout() {}
    
    /// Non-Trivial Constructor
    /**
        \param a_boxes          A Proto::Box representing each level of the hierarchy (coarse -> fine)
        \param a_periodicity    Defines the periodicity of the coarsest level's boundaries.
        \param a_split          If true, each domain will be split and load balanced according to 
                                MAXBOXSIZE and BLOCKFACTOR.
    */
    AMRLayout(std::vector<Proto::Box>& a_boxes, Proto::Point a_periodicity, bool a_split = true)
    {
        define(a_boxes, a_periodicity, a_split);
    }
    /// Define
    /**
        Define or redefine an existing AMRLayout
        
        \param a_boxes          A Proto::Box representing each level of the hierarchy (coarse -> fine)
        \param a_periodicity    Defines the periodicity of the coarsest level's boundaries.
        \param a_split          If true, each domain will be split and load balanced according to 
                                MAXBOXSIZE and BLOCKFACTOR.

    */
    void define(std::vector<Proto::Box>& a_boxes, Proto::Point a_periodicity, bool a_split = true)
    {
        PR_TIME("AMRLayout::define");
        m_layouts.resize(a_boxes.size());
        Proto::Box B = a_boxes[0];
        for (int ii = 0; ii < a_boxes.size(); ii++)
        {
            Proto::Box b = a_boxes[ii];
            Proto::Point p = a_periodicity;
            if (ii > 0) {
                int v[DIM];
                Proto::Box bc = b.coarsen(ipow(AMR_REFRATIO,ii));
                for (int d = 0; d < DIM; d++)
                {
                    if ((bc.low()[d] == B.low()[d]) && (bc.high()[d] == B.high()[d]))
                    {
                        v[d] = 1;
                    } else {
                        v[d] = 0;
                    }
                }
                p = Proto::Point(v);
            }
            buildLayout(m_layouts[ii], b, p, a_split);
        }
    }
    
    /// Level Layout Accessor
    /**
        Retrieve a reference to the Chombo::DisjointBoxLayout associated with level a_level.

        \param a_level  A level in the hierarchy (level 0 is coarsest)
    */
    const DisjointBoxLayout& operator[](int a_level) const {
        return m_layouts[a_level];
    }
   
    /// Level Domain Accessor
    /**
        Retrieve a reference to the Proto::Box representing the domain of level a_level.
        This function will likely return a std::vector of Proto::Box once branching is implemented.

        \param a_level  A level in the hierarchy (level 0 is coarsest)
    */
    Proto::Box domain(int a_level) const {
        return m_layouts[a_level].physDomain().domainBox();
    }
  
    /// Periodic Query
    /**
        Returns true if the layout is periodic in direction a_dir on level a_level.

        \param a_level  A level in the hierarchy (level 0 is coarsest)
        \param a_dir    An integer direction in [0,DIM)
    */ 
    bool periodic(int a_level, unsigned int a_dir) const {
        PROTO_ASSERT(a_dir < DIM, "a_dir should be less than DIM (= %i)", DIM);
        return m_layouts[a_level].physDomain().isPeriodic(a_dir);
    }
   
    /// Number of Levels
    /**
        Returns the number of levels in this layout.
    */
    int size() const {return m_layouts.size();}
    
private:
    std::vector<DisjointBoxLayout> m_layouts; ///< Chombo::DisjointBoxLayouts of the hierarchy
};

