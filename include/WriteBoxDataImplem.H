#include <cstdio>
#include <iostream>


using namespace std;

inline const char* WriteBoxData(BoxData<double>& a_array,
                    double dx)
{
    static int fileCount = 0;
    static char nameBuffer[10];

    if (!a_array.defined())
    {
        return nameBuffer;
    }

    sprintf(nameBuffer, "md%06d",fileCount);
    WriteBoxData(nameBuffer, a_array, dx);

    fileCount++;

    return nameBuffer;
};

inline void WriteBoxData(const char*           a_filename,
             BoxData<double>&  a_array,
             double dx)
{
    if (!a_array.defined()){return;}
    const char* const vars[] = {a_filename};
    double origin[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        origin[ii] = 0.0;
    }
  
    WriteBoxData(a_filename, a_array, vars, origin, dx);
}
/*
template<class T, unsigned int C, unsigned char D, unsigned char E>
void WriteBoxData(const char*                a_filename,
             BoxData<T, C, D, E>&   a_array,
             const char* const          a_varnames[C*D*E],
             const double*              a_corner,
             const double&              a_h)
{
    if (!a_array.defined()){return;}
    if (a_filename == NULL)
    {
        return;
    }
    //FILE* fp = vtk_open_file(a_filename);
    WriteBoxData(a_filename, a_array, a_varnames, a_corner, a_h);
    //vtk_close_file(fp);
}
*/

template<class T, unsigned int C, unsigned char D, unsigned char E>
inline void WriteBoxData(const char*       a_filename,
             BoxData<T,C,D,E>&      a_array,
             const char* const      a_varnames[C*D*E],
             const double*          a_corner,
             const double&          a_h)
{


    int dim[3] = {1,1,1};

    int vardims[C*D*E];
    int centering[C*D*E];
    double* vars[C*D*E];
    double* x[3] {nullptr,nullptr,nullptr};

    Point lo, hi;
    const Bx& box = a_array.box();
    lo = box.low();
    hi = box.high();
    for (int i = 0; i < DIM; i++)
    {
        dim[i] = hi[i] - lo[i] + 2;
        x[i] = (double*)StackAlloc::alloc((dim[i]+1)*sizeof(double));
        for (int k = 0; k < dim[i]+1; k++)
        {
            x[i][k] = (double)(a_corner[i] + a_h * k);
        }
    }

    if (DIM == 2)
    {
        x[2] = (double*)StackAlloc::alloc(2*sizeof(double));
        x[2][0] = 0.0;
        x[2][1] = 1.0;
    }

    for (int i = 0; i < C*D*E; i++)
    {
        vardims  [i] = 1;
        centering[i] = 0;

//        vars[i] = &((a_array(a_array.box().low(),i));
        vars[i] = a_array.data(a_array.box().low(),i);
    }

    write_rectilinear_mesh( a_filename,
                            dim,
                            x[0],x[1],x[2],
                            C*D*E,
                            vardims,
                            centering,
                            a_varnames,
                            vars);

    StackAlloc::free(x[2]); // it's stack, have to do it backwards :-)
    StackAlloc::free(x[1]);
    StackAlloc::free(x[0]);
}


