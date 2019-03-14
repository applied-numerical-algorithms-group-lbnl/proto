#pragma once

#define AMR_REFRATIO 2
#define MG_REFRATIO 2

#define PRE_RELAX 4
#define POST_RELAX 4
#define BOTTOM_RELAX 20

#define MAXBOXSIZE 8
#define BLOCKFACTOR 2

#include "AMRIO.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DisjointBoxLayout.H"

//using namespace Proto;

int AMR_LEVEL = 0;
int MG_LEVEL = 0;

template<typename T>
void printLevel(const LevelData<BaseFab<T>>& a_data)
{
    auto iter = a_data.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        const Proto::BoxData<T>& bd = a_data[iter];
        bd.printData();
    }
}

void printLevel(const LevelData<FArrayBox>& a_data)
{
    auto iter = a_data.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        const Proto::BoxData<Real>& bd = a_data[iter];
        bd.printData();
    }
}

void printFirst(const LevelData<FArrayBox>& a_data, int a_prec = 2)
{
    auto iter = a_data.dataIterator();
    iter.begin();
    const Proto::BoxData<Real>& patch = a_data[iter];
    patch.printData(a_prec);
}

void printBoxes(const DisjointBoxLayout& a_layout)
{
    auto iter = a_layout.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        std::cout << Proto::Box(a_layout[iter]) << std::endl;
    }
}

void buildLayout(DisjointBoxLayout& a_layout,
                 const Box& a_box,
                 const Proto::Point& a_point = Proto::Point::Ones(),
                 bool doSplit = true)
{
    bool periodic[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        periodic[ii] = (bool)a_point[ii];
    }
    ProblemDomain domain(a_box, periodic);
    Vector<Box> boxVect;
    Vector<int> procAssign;
    if (doSplit)
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

Real absMax(std::vector<std::shared_ptr<LevelData<FArrayBox>>>& a_data,
         bool a_interior = true,
         int a_edgeSize = 0)
{
    Real maxVal = 0;
    
    LevelData<FArrayBox>& topData = *(a_data[0]);
    auto topLayout = topData.disjointBoxLayout();
    Proto::Box topDomain = topLayout.physDomain().domainBox();

    Proto::Box topInterior = topDomain.grow(-a_edgeSize);
    auto titer = topLayout.dataIterator();
    for (titer.begin(); titer.ok(); ++titer)
    {
        Proto::BoxData<Real> patch = topData[titer];
        Proto::Box patchBox = topLayout[titer];
        for (auto iter = patchBox.begin(); iter != patchBox.end(); ++iter)
        {
            if ((a_interior && topInterior.contains(*iter))
                    || (!(a_interior) && !topInterior.contains(*iter)))
                maxVal = max(maxVal, abs(patch(*iter)));
        }
    }
    
    for (int ii = 1; ii < a_data.size(); ii++)
    {
        LevelData<FArrayBox>& fineData = *(a_data[ii]);
        LevelData<FArrayBox>& coarseData = *(a_data[ii-1]);
        auto coarseLayout = coarseData.disjointBoxLayout();
        auto fineLayout = fineData.disjointBoxLayout();
        Proto::Box coarseDomain = coarseLayout.physDomain().domainBox();
        Proto::Box coarseInterior = coarseDomain.grow(-a_edgeSize);
        Proto::Box fineDomain = fineLayout.physDomain().domainBox();
        Proto::Box coarseFineDomain = fineDomain.coarsen(AMR_REFRATIO);
        Proto::Box coarseFineBoundary = coarseFineDomain.grow(a_edgeSize);
        
        auto citer = coarseData.dataIterator();
        for (citer.begin(); citer.ok(); ++citer)
        {
            Proto::BoxData<Real> patch = coarseData[citer];
            Proto::Box patchBox = coarseLayout[citer];
            for (auto iter = patchBox.begin(); iter != patchBox.end(); ++iter)
            {
                if (!coarseFineDomain.contains(*iter))
                {
                    if ((a_interior && !coarseFineBoundary.contains(*iter) && coarseInterior.contains(*iter))
                            || (!(a_interior) && coarseFineBoundary.contains(*iter))
                            || (!(a_interior) && coarseDomain.contains(*iter) && !coarseInterior.contains(*iter)))
                    {
                        maxVal = max(maxVal, abs(patch(*iter)));
                    }
                }
            }
        }
    }
    return maxVal;
}


Real integrate(const LevelData<FArrayBox>& a_data, Real dx)
{
    auto iter = a_data.dataIterator();
    auto layout = a_data.disjointBoxLayout();
    Real sum = 0.0;
    for (iter.begin(); iter.ok(); ++iter)
    {
        Proto::BoxData<Real> patch = a_data[iter];
        Proto::Box b = layout[iter];
        for (auto biter = b.begin(); biter != b.end(); ++biter)
        {
            sum += patch(*biter)*pow(dx,DIM);
        }
    }
    return sum;
}

Real integrate(const std::vector<std::shared_ptr<LevelData<FArrayBox>>>& a_data, Real dx)
{
    Real sum = 0;
    for (int ii = 1; ii < a_data.size(); ii++)
    {
        LevelData<FArrayBox>& fineData = *(a_data[ii]);
        LevelData<FArrayBox>& coarseData = *(a_data[ii-1]);
        auto coarseLayout = coarseData.disjointBoxLayout();
        auto fineLayout = fineData.disjointBoxLayout();
        Proto::Box coarseDomain = coarseLayout.physDomain().domainBox();
        Proto::Box fineDomain = fineLayout.physDomain().domainBox();
        Proto::Box coarseFineDomain = fineDomain.coarsen(AMR_REFRATIO);
        
        auto citer = coarseData.dataIterator();
        for (citer.begin(); citer.ok(); ++citer)
        {
            Proto::BoxData<Real> patch = coarseData[citer];
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
                Proto::BoxData<Real> patch = fineData[fiter];
                Proto::Box patchBox = fineLayout[fiter];
                for (auto iter = patchBox.begin(); iter != patchBox.end(); ++iter)
                {
                    sum += patch(*iter)*dx*dx;
                }
            }
        }
    }
    return sum;

}

template<typename... Args>
void writeLevel(const LevelData<FArrayBox>& a_data, char* a_fname, Args... a_params)
{
    char fname[100];
    sprintf(fname, a_fname, a_params...);
    writeLevelname(&a_data, fname);

}

