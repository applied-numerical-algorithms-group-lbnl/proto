#pragma once

#include "Proto.H"
#include "ProtoInterface.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "AMRIO.H"
#include "AMRLayout.H"

/// Level AbsMax
/**
    Computes the maximum absolute value on a level

    \param a_data   A Chombo LevelData<FArrayBox> reference
*/
Real absMax(const LevelData<FArrayBox>& a_data)
{
    PR_TIME("::absMax(Chombo::LevelData) (utility found in AMRData.H)");
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
    PR_TIME("::integrate(Proto::BoxData) (utility found in AMRData.H)");
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
    PR_TIME("::integrate(Chombo::LevelData) (utility found in AMRData.H)");
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

/// Adaptive Mesh Refinement Hierarchy Data
/**
    This is the principle data holder object for the AMRFAS framework.
    
    Implementation details:
        At the moment, each level is implemented using a Chombo::LevelData<Chombo::FArrayBox>
        This implementation is well hidden from the user, so moving towards something like
        Proto::LevelData<Proto::BoxData<...>> should not affect user codes.
        Should this particular move take place, AMRData will use the same template parameters as
        Proto::BoxData
*/

template<int NUMCOMPS=1>
class AMRData {
public:
    
    /// Default Constructor
    AMRData(){m_layout = nullptr;}

    /// Non-Trivial Constructor
    /**
        This is the main constructor for AMRData and is recommended for public use

        \param a_layout     The AMRLayout upon which this is based
        \param a_ghost      Defines the number of ghost cells in each direction
        \param a_dx         Grid spacing (isotropic for now)
        \param a_initZero   If true, this will be zero-initialized everywhere

    */
    AMRData(AMRLayout& a_layout, Proto::Point a_ghost, Real a_dx, bool a_initZero = true)
    {
        define(a_layout, a_ghost, a_dx, a_initZero);
    }
   
    /// Define
    /**
        This is the main version of define and is recommended for public use. 
        
        \param a_layout     The AMRLayout upon which this is based
        \param a_ghost      Defines the number of ghost cells in each direction
        \param a_dx         Grid spacing (isotropic for now)
        \param a_initZero   If true, this will be zero-initialized everywhere
    */ 
    void define(AMRLayout& a_layout, Proto::Point a_ghost, Real a_dx, bool a_initZero = true)
    {
        PR_TIME("AMRData::define");
        m_ghost = a_ghost;
        m_layout = &a_layout;
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
    
    //TODO: This constructor should probably be removed since it's semantics are not intuitive -CLG
    /// Non-Copying Clone Constructor
    /**
        Creates an UNINITIALIZED AMRData based on the parameters of a_data. The resulting AMRData
        will have the same layout, ghost cell allowance, and grid spacing, but NO DATA IS COPIED. 

        Required by the IMEXIntegrator API (for now)
    */
    AMRData(AMRData<NUMCOMPS>& a_data)
    {
        define(a_data);
    }
   
    //TODO: This function should probably be removed since it's semantics are not intuitive -CLG
    /// Non-Copying Clone Define
    /**
        Defines an UNINITIALIZED AMRData based on the parameters of a_data. The resulting AMRData
        will have the same layout, ghost cell allowance, and grid spacing, but NO DATA IS COPIED. 

        Not recommended for public use.

        Required by the IMEXIntegrator API (for now)

        \param a_data   An AMRData
    */
    void define(AMRData<NUMCOMPS>& a_data)
    {
        define(a_data.layout(), a_data.m_ghost, a_data.m_dx, false);
    }
    
    /// Alternative Define
    /**
        A version of define with an (unused) ncomps parameter. Necessary for compliance with the
        IMEXIntegrator API (for now). For internal use only.

        \param a_layout     An AMR Hierarchy Layout
        \param a_ncomps     Number of components (this is completely unused as NCOMPS is a template parameter)
        \param a_ghost      Defines the number of ghost cells in each direction
    */
    void define(AMRLayout& a_layout, int a_ncomps, Proto::Point a_ghost)
    {
        PROTO_ASSERT(NUMCOMPS == a_ncomps, "number of components doesn't match template parameter");
        define(a_layout, a_ghost, 1, false);
    }
   
    /// Set Grid Spacing
    /**
        Change the grid spacing of this. No data is altered.
        Part of the IMEXIntegrator interface. 
        
        \param a_dx New Grid spacing.
    */ 
    void setGridSpacing(Real a_dx){m_dx = a_dx;}

    /// Clear Data
    /**
        Zeros out all the data in this.
        Part of the IMEXIntegrator interface. 
    */
    void clear() {
        initialize(
            [=] PROTO_LAMBDA (Proto::Point& a_pt, Proto::Var<Real, NUMCOMPS>& a_data, Real a_dx)
            {
                for (int cc = 0; cc < NUMCOMPS; cc++)
                {
                    a_data(cc) = 0.0;
                }
            }, m_dx
        );
    }
    
    /// Get AMRLayout
    /**
        Get the AMRLayout of this. 
        Part of the IMEXIntegrator interface. 
    */
    AMRLayout& layout() { return *m_layout;}
    
    /// Get Number of Components
    /**
        Get NUMCOMPS .
        Part of the IMEXIntegrator interface. 
    */
    constexpr int ncomps() { return NUMCOMPS;}
    
    /// Get Number of Ghost Cells
    /**
        Part of the IMEXIntegrator interface. 
    */
    Proto::Point ghost() {return m_ghost;}

    /// Initialize Hierarchy
    /**
        Initialize this using Proto::forall style syntax. 

        \param a_F  Function with the signature: void a_F(Proto::Point& a_p, Proto::Var<Real, NUMCOMPS>& a_data, Real a_dx);
        \param a_dx Coarsest grid spacing in the hierarchy TODO: I don't remember why this parameter is necessary...
    */
    template<typename Func>
    void initialize(const Func& a_F, double a_dx)
    {
        PR_TIME("AMRData::initialize");
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
    
    /// Get Level
    /**
        Retrieve a reference to the Chombo::LevelData representing level a_level.
        This is mostly intended for internal use, but is clearly useful for general use.
        Users should be aware that the implementation Chombo::LevelData<FArrayBox> is very likely to change.

        \param a_level  Desired level in [0,this.size()) (The coarsest level is 0)
    */
    LevelData<FArrayBox>& operator[](int a_level)
    {
        return *(m_data[a_level]);
    }
    
    /// Get Level (Const)
    /**
        Retrieve a reference to the Chombo::LevelData representing level a_level.
        This is mostly intended for internal use, but is clearly useful for general use.
        Users should be aware that the implementation Chombo::LevelData<FArrayBox> is very likely to change.

        \param a_level  Desired level in [0,this.size()) (The coarsest level is 0)
    */
    const LevelData<FArrayBox>& operator[](int a_level) const
    {
        return *(m_data[a_level]);
    }

    /// Number of Levels
    /**
        Returns the number of levels in this layout.
    */
    int size() const {return m_data.size();}
  
    
    /// Data Copy on Intersection
    /**
        Copy the data of this into a_data on their intersection. Ghost cells of a_data are updated.

        Implementation Details:
            This function is just a wrapper for Chombo::LevelData::copyTo which is called on each level,
            hence it has all the associated bells and whistles. 

        \param a_data   Another AMRData
    */
    void copyTo(AMRData<NUMCOMPS>& a_data) const
    {
        for (int ii = 0; ii < size(); ii++)
        {
            (*(m_data[ii])).copyTo(a_data[ii]);
        }
    }
   
    /// Data Copy and Coarsen on Intersection
    /**
        Coarsen this and copy into a_data. The coarsening ratio is AMR_REFRATIO by default.

        Implementation Details:
            The copying is done by Chombo::LevelData::copyTo and hence has all the associated attributes.
        
        \param a_data   Another AMRData
    */
    void coarsenTo(AMRData<NUMCOMPS>& a_data, int a_ratio = 0) const
    {
        PR_TIME("AMRData::coarsenTo");
        if (a_ratio == 0){a_ratio = AMR_REFRATIO;}
        auto S_avg = Proto::Stencil<Real>::AvgDown(a_ratio);
        for (int ii = 0; ii < size(); ii++)
        {
            DisjointBoxLayout tempLayout;
            coarsen_dbl(tempLayout, (*m_layout)[ii], a_ratio);
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
   
    /// Compute Cell Average
    /**
        Initializes a_data with the 2nd order cell average of the data in this. 
        
        \param a_data   Another AMRData
    */
    void toCellAverage(AMRData<NUMCOMPS>& a_data) const
    {
        PR_TIME("AMRData::toCellAverage");
        auto L = Proto::Stencil<Real>::Laplacian();
        for (int ii = 0; ii < size(); ii++)
        {
            auto iter = (*m_layout)[ii].dataIterator();
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
    
    /// Increment
    /**
        Adds the contents of a_data - optionally scaled - to the data in this.
        Both AMRData objects must have compatable AMRLayouts. 

        Required by IMEXIntegrator interface.

        \param a_data   another AMRData
        \param a_scale  a scaling
    */
    void increment(const AMRData<NUMCOMPS>& a_data, Real a_scale = 1.0) 
    {
        PR_TIME("AMRData::increment");
        for (int ii = 0; ii < size(); ii++)
        {
            auto& thisLevel = *(m_data[ii]);
            auto& otherLevel = a_data[ii];
            auto iter = (*m_layout)[ii].dataIterator();
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

    /// Multiply by Scalar
    /**
        Multiplies all data in this by a scalar.

        \param a_scale A scale
    */

    void times(Real a_scale)
    {
        PR_TIME("AMRData::times");
        for (int ii = 0; ii < size(); ii++)
        {
            auto& thisLevel = *(m_data[ii]);
            auto iter = (*m_layout)[ii].dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> ths;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(ths, thisLevel[iter]);
                ths *= a_scale;
            }
        }
    }

    /// Write to HDF5
    /**
        Outputs the contents of this to a VisIt compatable HDF5 file. Uses printf-style variadic arguments.

        \param a_fname      A filename pattern
        \param a_params     Optional parameters to the filename pattern
    */

    template<typename... Args>
    inline void write(const char* a_fname, Args... a_params)
    {
        PR_TIME("AMRData::write");
        char fname[100];
        sprintf(fname, a_fname, a_params...);
        Vector<DisjointBoxLayout> grids(size());
        Vector<LevelData<FArrayBox>*> data(size());
        for (int ii = 0; ii < grids.size(); ii++)
        {
            grids[ii] = (*m_layout)[ii];
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
        Box domain = (*m_layout).domain(0);
        //Real dx = m_mg.op().dx()*ipow(AMR_REFRATIO, m_level);
        Real dt = 1.0;
        Real t = 0.0;
        Vector<int> ratios(size()-1, AMR_REFRATIO);
        int numLevels = size();
        WriteAMRHierarchyHDF5(fname, grids, data, varNames, domain, m_dx, dt, t, ratios, numLevels);
    }

    /// AbsMax
    /**
        Computes the maximum absolute value on a whole AMR Hierarchy. Ignores invalid cells. 
    */
    Real absMax()
    {
        PR_TIME("AMRData::absMax");
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
      Integrates over a whole AMR Hierarchy using piecewise constant quadrature. Ignores coarse cells covered by a finer level.  
    */
    Real integrate()
    {
        PR_TIME("AMRData::absMax");
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
   
    /// Sum
    /**
        Computes the sum over the entire hierarchy. Ignores coarse cells covered by a finer layer. 
    */ 
    Real sum()
    {
        PR_TIME("AMRData::absMax");
        Real sum = 0.0;
        int n = 0;
        for (int ii = 1; ii < size(); ii++)
        {
            auto& coarseLevel = *(m_data[ii-1]);
            auto fineDomain = (*m_layout).domain(ii);
            auto iter = coarseLevel.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                Proto::BoxData<Real, NUMCOMPS> patch;
                ProtoCh::aliasBoxData<Real, NUMCOMPS>(patch, coarseLevel[iter]);
                Proto::Box b = (*m_layout)[ii-1][iter];
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
    AMRLayout* m_layout;///< layout
    std::vector<std::shared_ptr<LevelData<FArrayBox>>> m_data;///< data holders
    Proto::Point m_ghost; ///< Ghost cell structure
    Real m_dx; ///< Coarsest grid spacing
};

/// Write Level
/**
    Utility function for writing a Chombo::LevelData to a VisIt compatable HDF5 file. 
    Uses printf-style variadic arguments.

    \param a_fname      A filename pattern
    \param a_params     Optional parameters to the filename pattern
*/
template<typename... Args>
void writeLevel(const LevelData<FArrayBox>& a_data, const char* a_fname, Args... a_params)
{
    PR_TIME("::writeLevel(Chombo::LevelData) (Utility found in AMRData.H)");
    char fname[100];
    sprintf(fname, a_fname, a_params...);
    LevelData<FArrayBox> temp(a_data.disjointBoxLayout(), 1, Proto::Point::Zeros());
    a_data.copyTo(temp);
    writeLevelname(&a_data, fname);
}
