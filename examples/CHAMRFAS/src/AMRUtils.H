#pragma once

#define AMR_REFRATIO 2
#define MG_REFRATIO 2

#define PRE_RELAX 4
#define POST_RELAX 4
#define BOTTOM_RELAX 20

#define MAXBOXSIZE 16
#define BLOCKFACTOR 2

#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DisjointBoxLayout.H"

/// Print LevelData
/**
    Prints a LevelData to the command line. Pretty useless for real data, but helpful for debugging.

    \param a_data   a Chombo LevelData<FArrayBox>
*/
void printLevel(const LevelData<FArrayBox>& a_data)
{
    auto iter = a_data.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        const Proto::BoxData<Real>& bd = a_data[iter];
        bd.printData();
    }
}

/// Print Boxes in Layout
/**
    Prints all boxes in a Chombo DisjointBoxLayout. Note that the layout knows nothing about ghost cells

    \param a_layout     a Chombo DisjointBoxLayout reference
*/

void printBoxes(const DisjointBoxLayout& a_layout)
{
    auto iter = a_layout.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        std::cout << "\t" << Proto::Box(a_layout[iter]) << std::endl;
    }
}

/// Build DisjointBoxLayout
/**
    Utility function for initializing a load balanced DisjointBoxLayout.
    
    \param a_layout     An uninitialized DisjointBoxLayout
    \param a_box        Domain box of the layout
    \param a_point      Defines the periodicity in each direction (0 for non-periodic, periodic otherwise)
    \param a_doSplit    (Optional: Default = true) If true, do domain splitting and load balancing
*/
void buildLayout(DisjointBoxLayout& a_layout,
                 const Box& a_box,
                 const Proto::Point& a_periodic = Proto::Point::Ones(),
                 bool  a_doSplit = true)
{
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
    ProblemDomain domain = a_fineLayout.physDomain();
    domain.coarsen(a_refRatio);

    Vector<Box> boxVect;
    domainSplit(domain, boxVect, MAXBOXSIZE, BLOCKFACTOR);
    Vector<int> procAssign(boxVect.size(), 0);
    LoadBalance(procAssign, boxVect);
    a_coarseLayout.define(boxVect, procAssign, domain);
}

/// Level AbsMax
/**
    Computes the maximum absolute value on a level

    \param a_data   A Chombo LevelData<FArrayBox> reference
*/
Real absMax(const LevelData<FArrayBox>& a_data)
{
    Real maxVal = 0;
    auto iter = a_data.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        Proto::BoxData<Real> patch = a_data[iter];
        maxVal = max(patch.absMax(), maxVal);
    }
    return maxVal;
}

/// AMR AbsMax
/**
    Computes the maximum absolute value on a whole AMR Hierarchy interior or coarse-fine boundary. Ignores invalid cells. 

    \param a_data   An AMR Hierarchy
    \param a_level  Starting level. Computation proceeds from here to the finest level available. Input 0 to parse the whole hierarchy.
    \param a_interior   If true, only parses interior cells. If false, only parses boundary cells. 
    \param a_edgeSize   Width of the coarse-fine boundary. If 0, the interior is the entire (valid) hierarchy. 
*/
Real absMax(std::vector<std::shared_ptr<LevelData<FArrayBox>>>& a_data,
        int a_level=0,
        bool a_interior = true,
        int a_edgeSize = 0)
{
    PROTO_ASSERT(a_level < a_data.size(), "missed the base case somehow");
    LevelData<FArrayBox>& data = *(a_data[a_level]);
    if (a_level == a_data.size() - 1) //finest level
    {
        Real maxValue = 0.0;
        if (a_interior)
        {
            auto iter = data.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real> patch = data[iter];
                maxValue = std::max(maxValue, patch.absMax());
            }
        }
        return maxValue;
    } else {
        Real maxValue = absMax(a_data, a_level+1, a_interior, a_edgeSize);
        LevelData<FArrayBox>& fineData = *(a_data[a_level+1]);
        Proto::Box invalid = fineData.disjointBoxLayout().physDomain().domainBox();
        invalid = invalid.coarsen(AMR_REFRATIO);
        Proto::Box bounds = invalid.grow(a_edgeSize);
       
        auto iter = data.dataIterator();
        auto layout = data.disjointBoxLayout();
        for (iter.begin(); iter.ok(); ++iter)
        {
            Proto::BoxData<Real> patch = data[iter];
            Proto::Box b = layout[iter];
            for (auto piter = b.begin(); piter != b.end(); ++piter)
            {
                if (a_interior && !bounds.contains(*piter))
                {
                    maxValue = std::max(maxValue, std::abs(patch(*piter)));
                } else if (!a_interior && bounds.contains(*piter) && !invalid.contains(*piter))
                {
                    maxValue = std::max(maxValue, std::abs(patch(*piter)));
                }
            }
        }
        return maxValue;
    }
}

/// Patch Integrate
/**
    Integrate over a Proto::BoxData<Real> using piecewise constant quadrature.

    \param a_data   A Proto::BoxData<Real> to integrate over
    \param a_dx     The grid spacing
*/

Real integrate(const Proto::BoxData<Real> a_data,
    Proto::Box domain, Real dx)
{
    Real sum = 0.0;
    for (auto biter = domain.begin(); biter != domain.end(); ++biter)
    {
        sum += a_data(*biter)*pow(dx,DIM);
    }
    return sum;
}

/// Level Integrate
/**
    Integrates over a LevelData<FArrayBox> using piecewise constant quadrature. 

    \param a_data   A Chombo LevelData<FArrayBox> to integrate over
    \param a_dx     The grid spacing
*/

Real integrate(const LevelData<FArrayBox>& a_data, Real dx)
{
    auto iter = a_data.dataIterator();
    auto layout = a_data.disjointBoxLayout();
    Real sum = 0.0;
    for (iter.begin(); iter.ok(); ++iter)
    {
        const Proto::BoxData<Real> patch = a_data[iter];
        Proto::Box b = layout[iter];
        for (auto biter = b.begin(); biter != b.end(); ++biter)
        {
            sum += patch(*biter)*pow(dx,DIM);
        }
    }
    return sum;
}

/// AMR Integrate
/**
    Integrates over a whole AMR Hierarchy using piecewise constant quadrature. Ignores invalid cells. 

    \param a_data   An AMR Hierarchy
    \param a_dx     The coarsest grid spacing in the hierarchy
*/
Real integrate(const std::vector<std::shared_ptr<LevelData<FArrayBox>>>& a_data, Real dx)
{
    Real sum = 0;
    for (int ii = 1; ii < a_data.size(); ii++)
    {
        const LevelData<FArrayBox>& fineData = *(a_data[ii]);
        const LevelData<FArrayBox>& coarseData = *(a_data[ii-1]);
        auto coarseLayout = coarseData.disjointBoxLayout();
        auto fineLayout = fineData.disjointBoxLayout();
        Proto::Box coarseDomain = coarseLayout.physDomain().domainBox();
        Proto::Box fineDomain = fineLayout.physDomain().domainBox();
        Proto::Box coarseFineDomain = fineDomain.coarsen(AMR_REFRATIO);
        
        auto citer = coarseData.dataIterator();
        for (citer.begin(); citer.ok(); ++citer)
        {
            const Proto::BoxData<Real> patch = coarseData[citer];
            Proto::Box patchBox = coarseLayout[citer];
            for (auto iter = patchBox.begin(); iter != patchBox.end(); ++iter)
            {
                if (!coarseFineDomain.contains(*iter))
                {
                    sum += patch(*iter)*pow(dx,DIM);
                }
            }
        }
        dx /= AMR_REFRATIO;

        if (a_data.size() - ii == 1)
        {
            auto fiter = fineLayout.dataIterator();
            for (fiter.begin(); fiter.ok(); ++fiter)
            {
                const Proto::BoxData<Real> patch = fineData[fiter];
                Proto::Box patchBox = fineLayout[fiter];
                sum += integrate(patch, patchBox, dx);
                /*
                for (auto iter = patchBox.begin(); iter != patchBox.end(); ++iter)
                {
                    sum += patch(*iter)*pow(dx,DIM);
                }
                */
            }
        }
    }
    return sum;

}

template<typename... Args>
void writeLevel(const LevelData<FArrayBox>& a_data, const char* a_fname, Args... a_params)
{
    char fname[100];
    sprintf(fname, a_fname, a_params...);
    LevelData<FArrayBox> temp(a_data.disjointBoxLayout(), 1, Proto::Point::Zeros());
    a_data.copyTo(temp);
    writeLevelname(&a_data, fname);
}
