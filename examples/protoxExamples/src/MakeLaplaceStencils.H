#ifndef _MAKE_LAPLACE_STENCILS_H
#define _MAKE_LAPLACE_STENCILS_H
int getIndex(Point a_pt)
{
  Point abspt;
  for (int dir = 0 ; dir < DIM ; dir++)
    {
      abspt[dir] = abs(a_pt[dir]);
    }
  if (abspt[1] > abspt[0])
    {
      int tmp = abspt[0];
      abspt[0] = abspt[1];
      abspt[1] = tmp;
    }
#if DIM==3            
  if (abspt[2] > abspt[1])
    {
      int tmp = abspt[1];
      abspt[1] = abspt[2];
      abspt[2] = tmp;
    }
  if (abspt[1] > abspt[0])
    {
      int tmp = abspt[0];
      abspt[0] = abspt[1];
      abspt[1] = tmp;
    }
#endif              
  int ind;
  
  if (abspt[0] == 0)
    {
      ind = 0;
    }
  else if (abspt[0] == 1)
    {
      ind = abspt[0] + abspt[1] + abspt[2];
    }
  else if (abspt[0] == 2)
    {
      if (abspt[1] == 0)
        {
          ind = 4;
        }
      else if (abspt[1] == 1)
        {
          ind = 4 + abspt[1] + abspt[2];
        }
      else if (abspt[1] == 2)
        {
          ind = 4+abspt[1] + 1 + abspt[2];
        }
    }
  return ind;
}
template <typename T>
void LaplaceStencils3D(Stencil<T>& a_sten,
                       T& a_diagCoef,
                       int& a_diamStencil,
                       int a_solverType)
{
  Stencil<double> sten,sten0;
  Array<double,10> coefs;
  if ((a_solverType < 4) && (a_solverType > 0))
    {
      if (a_solverType==1)
        {
          Box UnitBox(Point::Ones(-1),Point::Ones(1)); 
          Array<double,10> coef_0 = {-4.0,1.0/3.0,1.0/6.0,0,0,0,0,0,0};
          coefs = coef_0;
          a_diamStencil = 1;
        }

      else if (a_solverType==2)
        {
          Array<double,10> coef_1 {-64.0/15,7.0/15,1.0/10.0,1.0/30,0,0,0,0,0,0};
          coefs = coef_1;
          a_diamStencil = 1;
        }
      else if (a_solverType == 4)
        {
          Array<double,10> coef_2 = {-6848.0/1575,  //000
                                     776.0/1575,  //100
                                     256.0/2835,  //110
                                     64.0/1575,   //111
                                     -11.0/2100,  //200
                                     4.0/2025,    //210
                                     -16.0/14175, //211
                                     -11.0/16200, //220
                                     2.0/14175,  //221
                                     -1.0/7560  };//222
          coefs = coef_2;
          a_diamStencil = 2;
        }
      double stentot = 0;
      Box unitBox(Point::Ones(-a_diamStencil),Point::Ones(a_diamStencil));
      for (auto bit : unitBox)
        {
          Point pt = bit;
          int ind = getIndex(pt);
          if (coefs[ind] != 0)
            {
              sten = sten + coefs[ind]*Shift(pt);
              stentot += coefs[ind];
            }
        }
      a_sten = sten;
      a_diagCoef = -coefs[0];
    }
  if (a_solverType == 3)
    {
      Stencil<T> stenTemp;
      for (int dir = 0; dir < DIM; dir++)
        {
          Stencil<double> sten1d =
            (-1.0/12)*Shift(Point::Basis(dir,2))+
            (-1.0/12)*Shift(Point::Basis(dir,-2))+
            (4.0/3)*Shift(Point::Basis(dir,1))+
            (4.0/3)*Shift(Point::Basis(dir,-1))+
            (-5.0/2)*Shift(Point::Zeros());
          stenTemp+= sten1d;
        }
      a_sten = stenTemp; 
      a_diamStencil = 2;
      a_diagCoef = 2.5*DIM;
    }
  if (a_solverType == 0)
    {
      a_sten = Stencil<T>::Laplacian();
      a_diamStencil = 1;
      a_diagCoef = 2*DIM;
    }
  if (a_solverType > 4)
    {
      cout << "no 3D solvertype = " << a_solverType << endl;
      abort();
    }
}
template <typename T>
void LaplaceStencils2D(Stencil<T>& a_sten,
                       T& a_diagCoef,
                       int& a_diamStencil,
                       int a_solverType)
  
{
  if (a_solverType == 0)
    {
      a_sten =
        1.0*Shift(Point::Basis(0,1)) + 1.0*Shift(Point::Basis(0,-1)) +
        1.0*Shift(Point::Basis(1,1)) + 1.0*Shift(Point::Basis(1,-1)) +
        (-4.0)*Shift(Point::Zeros());
      a_diagCoef = 4.0;
      a_diamStencil = 1;
    }
  if (a_solverType == 1)
    {
      a_sten =
        (4./6.)*Shift(Point::Basis(0,1)) + (4./6.)*Shift(Point::Basis(0,-1)) +
        (4./6.)*Shift(Point::Basis(1,1)) + (4./6.)*Shift(Point::Basis(1,-1)) +
        (1./6.)*Shift(Point::Ones(1)) + (1./6.)*Shift(Point::Ones(-1)) +
        (1./6.)*Shift(Point::Basis(0,1) + Point::Basis(1,-1))+ 
        (1./6.)*Shift(Point::Basis(0,-1) + Point::Basis(1,1)) + 
        (-20.0/6.0)*Shift(Point::Zeros());
      a_diagCoef = (20.0/6.0);
      a_diamStencil = 1;
    }
  if (a_solverType == 2)
    {
      Stencil<T> stenTemp;
      for (int dir = 0; dir < DIM; dir++)
        {
          Stencil<double> sten1d =
            (-1.0/12)*Shift(Point::Basis(dir,2))+
            (-1.0/12)*Shift(Point::Basis(dir,-2))+
            (4.0/3)*Shift(Point::Basis(dir,1))+
            (4.0/3)*Shift(Point::Basis(dir,-1))+
            (-5.0/2)*Shift(Point::Zeros());
          stenTemp+= sten1d;
        }
      a_sten = stenTemp; 
      a_diamStencil = 2;
      a_diagCoef = 2.5*DIM;
    }
  if (a_solverType > 2)
    {
      cout << "no 2D solvertype = " << a_solverType << endl;
      abort();
    }
}                          
template <typename T>
Stencil<T> LaplaceStencil(int a_solverType)
{
  Stencil<T> sten;
  T diagcoef;
  int diamStencil;
  
#if DIM == 3
  LaplaceStencils3D(sten,diagcoef,diamStencil,a_solverType);
#endif
#if DIM == 2
  LaplaceStencils2D(sten,diagcoef,diamStencil,a_solverType);
#endif
  return sten;
}
#endif
