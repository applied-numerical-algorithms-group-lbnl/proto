#pragma once

#define AMR_REFRATIO 2
#define MG_REFRATIO 2

#define PRE_RELAX 4
#define POST_RELAX 4
#define BOTTOM_RELAX 20

#define MAXBOXSIZE 8
#define BLOCKFACTOR 2

#include "Proto.H"
#include "ProtoInterface.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "AMRIO.H"

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
    double maxVal = 0;
    auto iter = a_data.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        Proto::BoxData<Real> patch;
        ProtoCh::aliasBoxData<Real, 1>(patch, a_data[iter]);
        maxVal = max(patch.absMax(), maxVal);
    }
#ifdef CH_MPI
    std::cout << "local max: " << maxVal << std::endl;
    double globalMax = 0;
    MPI_Reduce(&maxVal, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return globalMax;
#endif
    return maxVal;
}

/// Patch Integrate
/**
    Integrate over a Proto::BoxData<Real> using piecewise constant quadrature.

    \param a_data   A Proto::BoxData<Real> to integrate over
    \param a_dx     The grid spacing
*/
Real integrate(const Proto::BoxData<Real>& a_data,
    Proto::Box domain, Real dx)
{
    double sum = 0.0;
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
    double sum = 0.0;
    for (iter.begin(); iter.ok(); ++iter)
    {
        Proto::BoxData<Real> patch;
        ProtoCh::aliasBoxData<Real, 1>(patch, a_data[iter]);
        Proto::Box b = layout[iter];
        sum += integrate(patch, b, dx);
    }
#ifdef CH_MPI
    double globalSum = 0.0;
    MPI_Reduce(&sum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return globalSum;
#endif
    return sum;
}

class AMRLayout {
public:
    AMRLayout() {}
    
    AMRLayout(std::vector<Proto::Box>& a_boxes, Proto::Point a_periodicity, bool a_split = true)
    {
        define(a_boxes, a_periodicity, a_split);
    }

    void define(std::vector<Proto::Box>& a_boxes, Proto::Point a_periodicity, bool a_split = true)
    {
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

    const DisjointBoxLayout& operator[](int a_level) const {
        return m_layouts[a_level];
    }
   
    Proto::Box domain(int a_level) const {
        return m_layouts[a_level].physDomain().domainBox();
    }
   
    bool periodic(int a_level, unsigned int a_dir) const {
        PROTO_ASSERT(a_dir < DIM, "a_dir should be less than DIM (= %i)", DIM);
        return m_layouts[a_level].physDomain().isPeriodic(a_dir);
    }
   
    int size() const {return m_layouts.size();}
    
private:
    
    std::vector<DisjointBoxLayout> m_layouts;
};

template<int NUMCOMPS=1>
class AMRData {
public:
    AMRData(AMRLayout& a_layout, Proto::Point a_ghost, Real a_dx, bool a_initZero = true)
    : m_layout(a_layout)
    {
        define(a_layout, a_ghost, a_dx, a_initZero);
    }

    void define(AMRLayout& a_layout, Proto::Point a_ghost, Real a_dx, bool a_initZero = true)
    {
        m_layout = a_layout;
        m_dx = a_dx;
        m_data.resize(a_layout.size());
        for (int ii = 0; ii < a_layout.size(); ii++)
        {
            m_data[ii] = make_shared<LevelData<FArrayBox>>(a_layout[ii], NUMCOMPS, a_ghost);
        }
        if (a_initZero)
        {
            initialize(
                [=] PROTO_LAMBDA (Proto::Point& a_pt, Proto::Var<Real, NUMCOMPS>& a_data, Real a_dx)
                {
                    for (int cc = 0; cc < NUMCOMPS; cc++)
                    {
                        a_data(cc) = 0.0;
                    }
                }, a_dx);
        }
    }

    /// Initialize Hierarchy
    /**
        \param a_F  Function with the syntax: void a_F(Proto::Point& a_p, Proto::Var<Real, NUMCOMPS>& a_data, Real a_dx);
        \param a_dx Coarsest grid spacing in the hierarchy
    */
    template<typename Func>
    void initialize(const Func& a_F, double a_dx)
    {
        m_dx = a_dx;
        for (int ii = 0; ii < size(); ii++)
        {
            double dx = a_dx/pow(AMR_REFRATIO,ii);
            auto& LD = *(m_data[ii]);
            auto iter = LD.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> patch;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, LD[iter]);
                forallInPlace_p(a_F, patch, dx);
            }
            LD.exchange();
        }
    }

    LevelData<FArrayBox>& operator[](int a_level)
    {
        return *(m_data[a_level]);
    }
    
    const LevelData<FArrayBox>& operator[](int a_level) const
    {
        return *(m_data[a_level]);
    }

    int size() const {return m_data.size();}
   
    void copyTo(AMRData<NUMCOMPS>& a_data) const
    {
        for (int ii = 0; ii < size(); ii++)
        {
            (*(m_data[ii])).copyTo(a_data[ii]);
        }
    }
   
    void coarsenTo(AMRData<NUMCOMPS>& a_data, int a_ratio = 0) const
    {
        if (a_ratio == 0){a_ratio = AMR_REFRATIO;}
        auto S_avg = Proto::Stencil<Real>::AvgDown(a_ratio);
        for (int ii = 0; ii < size(); ii++)
        {
            DisjointBoxLayout tempLayout;
            coarsen_dbl(tempLayout, m_layout[ii], a_ratio);
            LevelData<FArrayBox> tempData(tempLayout, NUMCOMPS, Proto::Point::Zeros());
            
            auto iter = tempLayout.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> tmp;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(tmp, tempData[iter]);
                Proto::BoxData<Real, NUMCOMPS> src;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(src, (*(m_data[ii]))[iter]);
                tmp |= S_avg(src);
            }
            tempData.copyTo(a_data[ii]);
        }
    }
   
    void toCellAverage(AMRData<NUMCOMPS>& a_data) const
    {
        auto L = Proto::Stencil<Real>::Laplacian();
        for (int ii = 0; ii < size(); ii++)
        {
            auto iter = m_layout[ii].dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> src;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(src, (*(m_data[ii]))[iter]);
                Proto::BoxData<Real, NUMCOMPS> dst;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(dst, a_data[ii][iter]);
                PROTO_ASSERT(L.range(src.box()).contains(dst.box()), "AMRData::toCellAverage destination box is larger than Laplacian range");
                dst |= L(src, 1.0/24.0);
                dst += src;
            }
        }
    }
    
    void add(const AMRData<NUMCOMPS>& a_data, Real a_scale = 1.0) 
    {
        for (int ii = 0; ii < size(); ii++)
        {
            auto& thisLevel = *(m_data[ii]);
            auto& otherLevel = a_data[ii];
            auto iter = m_layout[ii].dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> ths;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(ths, thisLevel[iter]);
                Proto::BoxData<Real, NUMCOMPS> othr;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(othr, otherLevel[iter]);
                Proto::BoxData<Real, NUMCOMPS> tmp(othr.box());
                othr.copyTo(tmp);
                tmp *= a_scale;
                ths += tmp;
            }
        }
    }

    template<typename... Args>
    inline void write(const char* a_fname, Args... a_params)
    {
        char fname[100];
        sprintf(fname, a_fname, a_params...);
        Vector<DisjointBoxLayout> grids(size());
        Vector<LevelData<FArrayBox>*> data(size());
        for (int ii = 0; ii < grids.size(); ii++)
        {
            grids[ii] = m_layout[ii];
            auto& temp = *(m_data[ii]);
            data[ii] = &temp;
        }

        Vector<string> varNames(NUMCOMPS);
        for (int ii = 0; ii < varNames.size(); ii++)
        {
            char name[100];
            sprintf(name, "var_%i", ii);
            std::string s(name);
            varNames[ii] = s;
        }

        //Box domain = (*(a_data[0])).disjointBoxLayout().physDomain().domainBox();
        Box domain = m_layout.domain(0);
        //Real dx = m_mg.op().dx()*ipow(AMR_REFRATIO, m_level);
        Real dt = 1.0;
        Real t = 0.0;
        Vector<int> ratios(size()-1, AMR_REFRATIO);
        int numLevels = size();
        WriteAMRHierarchyHDF5(fname, grids, data, varNames, domain, m_dx, dt, t, ratios, numLevels);
    }

    /// AbsMax
    /**
      Computes the maximum absolute value on a whole AMR Hierarchy interior or coarse-fine boundary. Ignores invalid cells. 

      \param a_data   An AMR Hierarchy
      \param a_level  Starting level. Computation proceeds from here to the finest level available. Input 0 to parse the whole hierarchy.
      \param a_interior   If true, only parses interior cells. If false, only parses boundary cells. 
      \param a_edgeSize   Width of the coarse-fine boundary. If 0, the interior is the entire (valid) hierarchy. 
    */
    
    Real absMax()
    {
        double maxValue = 0.0;
        for (int level = 1; level < size(); level++)
        {
            LevelData<FArrayBox>& coarseData = *(m_data[level-1]);
            LevelData<FArrayBox>& fineData   = *(m_data[level]);
            auto& coarseLayout = coarseData.disjointBoxLayout();
            auto& fineLayout   = fineData.disjointBoxLayout();
            Proto::Box invalid = fineLayout.physDomain().domainBox();
            invalid = invalid.coarsen(AMR_REFRATIO);
            auto iter = coarseData.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> patch;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, coarseData[iter]);
                Proto::Box B = coarseLayout[iter];
                for (auto biter = B.begin(); biter != B.end(); ++biter)
                {
                    if (!invalid.contains(*biter))
                    {
                        maxValue = std::max(patch(*biter), maxValue);
                    }
                }
            }
        }
        LevelData<FArrayBox>& finestData = *(m_data[size()-1]);
        auto& finestLayout = finestData.disjointBoxLayout();
        auto iter = finestData.dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            //BoxData patch = finestData[iter];

            Proto::BoxData<Real> patch;
            ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, finestData[iter]);
            
            Proto::Box B = finestLayout[iter];
            for (auto biter = B.begin(); biter != B.end(); ++biter)
            {
                maxValue = std::max(patch(*biter), maxValue);
            }
        }
#ifdef CH_MPI
        double globalMax = 0.0;
        MPI_Allreduce(&maxValue, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        return globalMax;
#endif
        return maxValue;
    }

    /// Integrate
    /**
      Integrates over a whole AMR Hierarchy using piecewise constant quadrature. Ignores invalid cells. 
    */
    Real integrate()
    {
        double sum = 0.0;
        Real dx = m_dx;
        for (int ii = 1; ii < size(); ii++)
        {
            const LevelData<FArrayBox>& fineData = *(m_data[ii]);
            const LevelData<FArrayBox>& coarseData = *(m_data[ii-1]);
            auto coarseLayout = coarseData.disjointBoxLayout();
            auto fineLayout = fineData.disjointBoxLayout();
            Proto::Box coarseDomain = coarseLayout.physDomain().domainBox();
            Proto::Box fineDomain = fineLayout.physDomain().domainBox();
            Proto::Box coarseFineDomain = fineDomain.coarsen(AMR_REFRATIO);

            auto citer = coarseData.dataIterator();
            for (citer.begin(); citer.ok(); ++citer)
            {
                Proto::BoxData<Real> patch;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, coarseData[citer]);
                Proto::Box patchBox(coarseLayout[citer]);
                for (auto iter = patchBox.begin(); iter != patchBox.end(); ++iter)
                {
                    if (!coarseFineDomain.contains(*iter))
                    {
                        sum += patch(*iter)*pow(dx,DIM);
                    }
                }
            }
            dx /= AMR_REFRATIO;

            if (size() - ii == 1)
            {
                auto fiter = fineLayout.dataIterator();
                for (fiter.begin(); fiter.ok(); ++fiter)
                {
                    Proto::BoxData<Real> patch;
                    ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, fineData[fiter]);
                    Proto::Box patchBox = fineLayout[fiter];
                    sum += ::integrate(patch, patchBox, dx);
                }
            }
        }
#ifdef CH_MPI
        double globalSum = 0.0;
        MPI_Reduce(&sum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        return globalSum;
#endif
        return sum;
    }
    
    Real sum()
    {
        Real sum = 0.0;
        int n = 0;
        for (int ii = 1; ii < size(); ii++)
        {
            auto& coarseLevel = *(m_data[ii-1]);
            auto fineDomain = m_layout.domain(ii);
            auto iter = coarseLevel.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> patch;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, coarseLevel[iter]);
                Proto::Box b = m_layout[ii-1][iter];
                for (auto biter = b.begin(); biter != b.end(); ++biter)
                {
                    if (!fineDomain.contains(*biter)){sum += patch(*biter);}
                }
            }

            if (ii == (size() - 1))
            {
                auto& fineLevel = *(m_data[ii]);
                auto fiter = fineLevel.dataIterator();
                for (fiter.begin(); fiter.ok(); ++fiter)
                {
                    Proto::BoxData<Real, NUMCOMPS> patch;
                    ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, fineLevel[fiter]);
                    sum += patch.sum();
                }
            }
        }
#ifdef CH_MPI
        double globalSum = 0.0;
        MPI_Reduce(&sum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        return globalSum;
#endif
        return sum;
    }
    
private:
    AMRLayout& m_layout;
    std::vector<std::shared_ptr<LevelData<FArrayBox>>> m_data;
    Real m_dx; ///< Coarsest grid spacing
};

template<typename... Args>
void writeLevel(const LevelData<FArrayBox>& a_data, const char* a_fname, Args... a_params)
{
    char fname[100];
    sprintf(fname, a_fname, a_params...);
    LevelData<FArrayBox> temp(a_data.disjointBoxLayout(), 1, Proto::Point::Zeros());
    a_data.copyTo(temp);
    writeLevelname(&a_data, fname);
}
