#ifndef _ADVECTION_TEST_H_
#define _ADVECTION_TEST_H_
#include "Proto.H"
#include "AdvectionOp.H"
PROTO_KERNEL_START
void
initializeF(
        const Point& a_pt,
        Var<double,NUMCOMPS>& a_U,
        const double& a_h,
        const double& a_time)
{
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        a_U(comp) = 1.;
        for (int dir = 0; dir < DIM ; dir++)
        {
            a_U(comp)*=(-cos(2*M_PI*(a_pt[dir]*a_h - a_time) + 2*M_PI*a_h) + cos(2*M_PI*(a_pt[dir]*a_h - a_time)))/(2*M_PI*a_h);
        }
    }
}
PROTO_KERNEL_END(initializeF,initialize);
void
compactF(
         const Point& a_pt,
         Var<double,NUMCOMPS>& a_U,
         const double& a_h,
         const double& a_time)
{
  double r0 = .125;
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
      double rsq = 0.;
        for (int dir = 0; dir < DIM ; dir++)
        {
          double xcen = fmod(.5 + a_time,1.); 
          double xdir = a_pt[dir]*a_h + .5*a_h;
          double del; 
          if (xcen > xdir) del = min(abs(xdir - xcen),abs(xdir - xcen + 1.));
          if (xcen <= xdir) del = min(abs(xdir - xcen),abs(xdir - xcen - 1.));
          rsq += pow(del,2);
            //rsq += pow(a_pt[dir]*a_h + .5*a_h - .5 - a_time,2);
        }
        double r = sqrt(rsq);
        if (r > r0)
          {
            a_U(comp) = 1.;
          }
        else
          {
            a_U(comp) = 1.+pow(cos(M_PI*r/r0/2),6);
          }
    }
}
PROTO_KERNEL_END(compactF,compact);
template<typename T>
void advectionExact(
        LevelBoxData<T,NUMCOMPS>& a_data,
        T a_dx,
        T a_time)
{
    for (auto dit = a_data.begin();*dit != dit.end();++dit)
    {      
      //forallInPlace_p(initialize,a_data[*dit],a_dx,a_time);
      BoxData<T,NUMCOMPS> pointwise(a_data[*dit].box().grow(Point::Ones()));
      forallInPlace_p(compact,pointwise,a_dx,a_time);
      BoxData<T,NUMCOMPS> avCorr = Stencil<double>::Laplacian()(pointwise,1.0/24);
      pointwise.copyTo(a_data[*dit]);
      a_data[*dit] += avCorr;
    }
}

template<typename T>
void advectionError(
        LevelBoxData<T,NUMCOMPS>& a_error,
        const LevelBoxData<T,NUMCOMPS>& a_U,
        const T& a_dx,
        const T& a_dt,
        const T& a_time)
{

    advectionExact<T>(a_error,a_dx,a_time);
    for (auto dit = a_error.begin();*dit != dit.end();++dit)
    {
        a_error[*dit] -= a_U[*dit];
        a_error[*dit] *= (1./a_dt);
    }
}

template<typename T>
void advectionErrorPatch(
        BoxData<T,NUMCOMPS>& a_error,
        const BoxData<T,NUMCOMPS>& a_U,
        const T& a_dx,
        const T& a_dt,
        const T& a_time)
{

    forallInPlace_p(initialize,a_error,a_dx,a_time);
    a_error -= a_U;
    a_error *= (1./a_dt);
}
PROTO_KERNEL_START
    void
divergenceF(
        const Point& a_pt,
        Var<double,NUMCOMPS>& a_divF,
        const double& a_h,
        const double& a_time)
{
    // Notice that we are computing -div(F).
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        a_divF(comp) = 1.;
        for (int dir = 0; dir < DIM ; dir++)
        {
            a_divF(comp)*=-(sin(2*M_PI*(a_pt[dir]*a_h - a_time) +  2*M_PI*a_h)
                    - sin(2*M_PI*(a_pt[dir]*a_h - a_time)))/a_h;
        }

    }
#if 0
    if (a_pt[1] == 0)
    {
        cout << a_pt[0] << " , " << a_divF(0) << endl;
    }
#endif
}
PROTO_KERNEL_END(divergenceF,divergence);
    template<typename T>
void advectionDivergenceExact(
        LevelBoxData<T,NUMCOMPS>& a_data,
        const T& a_dx,
        const T& a_time)
{
    for (auto dit = a_data.begin();*dit != dit.end();++dit)
    {
        forallInPlace_p(divergence,a_data[*dit],a_dx,a_time);
    }
}
    template<typename T>
void writeDivergenceError(
        const LevelBoxData<T,NUMCOMPS>& a_data,
        const T& a_dx,
        const T& a_dt,
        const T& a_time,
        const string& a_str)
{
    LevelBoxData<T,NUMCOMPS> error(a_data.layout(),Point::Zeros());
    for (auto dit = a_data.begin();*dit != dit.end();++dit)
    {
        Box bx = a_data.layout()[*dit];
        if (a_dt > 0)
        {
            BoxData<T,NUMCOMPS> newtimedata(bx);
            T newtime = a_time+a_dt;
            forallInPlace_p(initialize,error[*dit],a_dx,a_time);
            forallInPlace_p(initialize,newtimedata,a_dx,newtime);
            error[*dit] -= newtimedata;
            error[*dit] += a_data[*dit];
            error[*dit] /= a_dt;
        }
        else
        {
            forallInPlace_p(divergence,error[*dit],a_dx,a_time);
            error[*dit]-=a_data[*dit];
        }    
    }
    //HDF5Handler h5;
    //h5.writeLevel(a_dx,error,a_str);
}
#endif
