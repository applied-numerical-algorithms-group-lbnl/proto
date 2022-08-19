//#include "Proto.H"
//#include "Proto_WriteBoxData.H"
#include "Proto_VisitWriter.H"
//#include "writers.h"
//#if 0
inline void write_point_mesh2(const char *filename, int npts, double *pts,
                             int nvars, int *vardim, const char *const *varnames,
                             double **vars)
{
  FILE *fp = vtk_open_file(filename);
  int numInColumn = 0;

  int i;
  char str[128];
  int *centering = NULL;

  write_header(fp);

  write_string("DATASET UNSTRUCTURED_GRID\n", fp);
  sprintf(str, "POINTS %d double\n", npts);
  write_string(str, fp);
  for (i = 0; i < 3 * npts; i++)
  {
    write_double(pts[i], fp, numInColumn);
  }
  return;
  new_section(fp, numInColumn);
  sprintf(str, "CELLS %d %d\n", npts, 2 * npts);
  write_string(str, fp);

  for (i = 0; i < npts; i++)
  {
    write_int(1, fp, numInColumn);
    write_int(i, fp, numInColumn);
    end_line(fp, numInColumn);
  }
  centering = (int *)malloc(nvars * sizeof(int));
  for (i = 0; i < nvars; i++)
    centering[i] = 1;
  write_variables(nvars, vardim, centering, varnames, vars, npts, npts, fp, numInColumn);
  free(centering);

  vtk_close_file(fp, numInColumn);
}
//#endif
template<unsigned int NPCOMP>
inline void PWrite(
            const vector<array<double,NPCOMP> >& a_data,
            const char *const *a_varNames,
            const string& a_fileName,
            int a_plotNumber)
{
  vector<vector<double>> vars(NPCOMP);
  unsigned int size = a_data.size();
  vector<double> x(3 * size);
  for (int icomp=0; icomp< NPCOMP; icomp++)
    {  
      vars[icomp] = vector<double>(size);
    }
  for (unsigned int i = 0; i < size; ++i)
  {
    for (int icomp = 0; icomp < NPCOMP;icomp++)
      {
        vars[icomp][i] = a_data[i][icomp];
      }
    x[i * 3] = a_data[i][0];
    x[i * 3 + 1] = a_data[i][1];
    x[i * 3 + 2] = 0.;
  }
  double *varPtr[NPCOMP];
  int vardim[NPCOMP];
  for (int icomp = 0; icomp < NPCOMP;icomp++)
    {
      varPtr[icomp] = &vars[icomp][0];
      vardim[icomp] = 1;
      
    }
  auto curFilename = a_fileName + to_string(a_plotNumber);
  int len = curFilename.length();
  char char_array[len + 1];
  strcpy(char_array, curFilename.c_str());
  write_point_mesh2(char_array, size, &(x[0]), NPCOMP, vardim, a_varNames, varPtr);
}
