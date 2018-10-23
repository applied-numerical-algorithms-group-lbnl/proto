#pragma once

#define AMR_REFRATIO 2
#define MG_REFRATIO 2

#define PRE_RELAX 4
#define POST_RELAX 4
#define BOTTOM_RELAX 20

#define MAXBOXSIZE 8
#define BLOCKFACTOR 2

#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "DisjointBoxLayout.H"

using namespace Proto;

void buildLayout(DisjointBoxLayout& a_layout,
                 const Bx& a_box,
                 const Point& a_point = Point::Ones()) //defaults to periodic
{
    bool periodic[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        periodic[ii] = (bool)a_point[ii];
    }
    ProblemDomain domain(a_box, periodic);
    Vector<Box> boxVect;
    domainSplit(domain, boxVect, MAXBOXSIZE, BLOCKFACTOR);

    Vector<int> procAssign(boxVect.size(), 0);

    LoadBalance(procAssign, boxVect);
    a_layout.define(boxVect, procAssign, domain);
}

void coarsenGrid(DisjointBoxLayout& a_coarseLayout,
                 const DisjointBoxLayout& a_fineLayout)
{
    ProblemDomain domain = a_fineLayout.physDomain();
    domain.coarsen(MG_REFRATIO);

    Vector<Box> boxVect;
    domainSplit(domain, boxVect, MAXBOXSIZE, BLOCKFACTOR);
    Vector<int> procAssign(boxVect.size(), 0);
    LoadBalance(procAssign, boxVect);
    a_coarseLayout.define(boxVect, procAssign, domain);
}
