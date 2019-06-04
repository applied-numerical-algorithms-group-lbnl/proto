#pragma once

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "AMRIO.H"

class AMRLayout {
public:
    AMRLayout(std::vector<Proto::Box>& a_boxes, Proto::Point a_periodicity, bool a_split = true)
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
                std::cout << "coarsest Box: " << B << ", this box: " << bc << std::endl;
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
            this->buildLayout(m_layouts[ii], b, p, a_split);
        }
    }

    DisjointBoxLayout& operator()(int a_level) {
        return m_layouts[a_level];
    }
   
    Proto::Box domain(int a_level) {
        return m_layouts[a_level].physDomain().domainBox();
    }
   
    bool periodic(int a_level, unsigned int a_dir) {
        PROTO_ASSERT(a_dir < DIM, "a_dir should be less than DIM (= %i)", DIM);
        return m_layouts[a_level].physDomain().isPeriodic(a_dir);
    }
   
    int size() {return m_layouts.size();}
    
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
private:
    
    std::vector<DisjointBoxLayout> m_layouts;
};

template<int NUMCOMPS, class DATA=FArrayBox>
class AMRHierarchy {
public:
    AMRHierarchy(AMRLayout& a_layout, Proto::Point a_ghost, Real a_dx, bool a_initZero = true)
    : m_layout(a_layout), m_dx(a_dx)
    {
        m_data.resize(a_layout.size());
        for (int ii = 0; ii < a_layout.size(); ii++)
        {
            m_data[ii] = make_shared<LevelData<DATA>>(a_layout(ii), NUMCOMPS, a_ghost);
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
                Proto::BoxData<Real, NUMCOMPS> patch = LD[iter];
                forallInPlace_p(a_F, patch, dx);
            }
        }
    }

    LevelData<DATA>& operator()(int a_level)
    {
        return *(m_data[a_level]);
    }

    int size() {return m_data.size();}
    
    template<typename... Args>
        inline void write(const char* a_fname, Args... a_params)
        {
            char fname[100];
            sprintf(fname, a_fname, a_params...);
            Vector<DisjointBoxLayout> grids(size());
            Vector<LevelData<FArrayBox>*> data(size());
            for (int ii = 0; ii < grids.size(); ii++)
            {
                grids[ii] = m_layout(ii);
                auto& temp = this->operator()(ii);
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

private:
    AMRLayout& m_layout;
    std::vector<std::shared_ptr<LevelData<DATA>>> m_data;
    Real m_dx; ///< Coarsest grid spacing
};

