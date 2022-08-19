#ifndef WRITERS_H
#define WRITERS_H

#include <vector>

using namespace std;

inline void write_point_mesh2(const char *, int, double *, int , int *, const char *const *, double **);

template<unsigned int NPCOMP>
inline void PWrite(
                   const vector<array<double,NPCOMP> >& a_data,
                   const char *const* a_varnames,
                   const string& a_filename,
                   int a_plotNumber);
#include "writersImplem.H"

#endif
