#ifndef _PROTO_WRITERECTMDARRAY_H_
#define _PROTO_WRITERECTMDARRAY_H_
#include "BoxData.H"
#include "VisitWriter.H"
#include "Stack.H"

#include <cstring>
#include <cstdio>

// first two are handy when debugging code.  Can call the first from debugger
using namespace Proto;

inline const char* WriteBoxData(BoxData<double>& a_array,
                    double dx = 1.0);


inline void WriteBoxData(const char*                a_filename,
                   BoxData<double>& a_array,
                   double dx = 1.0);


template<class T, unsigned int C, unsigned char D, unsigned char E>
inline void WriteBoxData(const char*           a_filename,
             BoxData<T, C, D, E>&       a_array,
             const char* const          a_varnames[C*D*E],
             const double*              a_corner,
             const double&              a_h);

template<class T, unsigned int C, unsigned char D, unsigned char E>
inline void WriteBoxData(FILE*                 a_fp,
             BoxData<T,C,D,E>&          a_array,
             const char* const          a_varnames[C*D*E],
             const double*              a_corner,
             const double&              a_h);

using std::string;
//convenience function, turn off iter in filename by setting a_iter < 0
template <int n>
inline void WriteData(BoxData<double,n> &    a_state,
               const   int       &    a_iter, 
               const   double    &    a_dx,
               const   string    &    a_varnamebase,
               const   string    &    a_filenamebase)
{
  double origin[DIM];
  for (int ii = 0; ii < DIM; ii++)
  {
    origin[ii] = 0.0;
  }

  char basename[1024];
  if(a_iter >= 0) 
  {
    sprintf(basename,"%s.%06d",a_filenamebase.c_str(),a_iter);
  }
  else
  {
    sprintf(basename,"%s",a_filenamebase.c_str());
  }
    
  char* velnames[n];
  for(int idir = 0; idir < n; idir++)
  {
    velnames[idir] = new char[1024];
    sprintf(velnames[idir],"%s%06d",a_varnamebase.c_str(),idir);
  }

  WriteBoxData(basename,a_state,velnames,origin,a_dx);

  for(int idir = 0; idir < n; idir++)
  {
    delete[] velnames[idir];
  }

};  

template <int n> inline void viewData(BoxData<double, n>* a_statePtr)
{
  WriteData<n>(*a_statePtr, -1, 1.0, string("var"), string("tempdata"));
  system("visit -o tempdata.vtk");
}

#include "WriteBoxDataImplem.H"

#endif 
