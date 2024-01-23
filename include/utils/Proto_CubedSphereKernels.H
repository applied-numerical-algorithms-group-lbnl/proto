#ifndef _PROTO_CUBED_SPHERE_KERNELS_
#define _PROTO_CUBED_SPHERE_KERNELS_
#define CUBED_SPHERE_SHELL_R0 1.0
#define CUBED_SPHERE_SHELL_R1 1.5
#define CUBED_SPHERE_SHELL_RADIAL_COORD 0
typedef Proto::Array<Proto::Array<double,DIM>,DIM> matType;
typedef Proto::Array<double,DIM> vecType;
using namespace Operator;
//==============================================================
// FORALL KERNELS

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_matMatProd_F(
                    Var<T,DIM,MEM,DIM>& a_prod,
                    const Var<T,DIM,MEM,DIM>& a_left,
                    const Var<T,DIM,MEM,DIM>& a_right)
{
  for (int ind0 = 0; ind0 < DIM; ind0++)
    {
      for (int ind1 = 0; ind0 < DIM; ind1++)
        {      
          a_prod(ind0,ind1) =0.;
          for (int ind = 0; ind < DIM; ind++)
            {
              a_prod(ind0,ind1) +=
                a_left(ind0,ind)*a_right(ind,ind1);
            }
        }
    }
}
PROTO_KERNEL_END(f_matMatProd_F,f_matMatProd)

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_matVecProd_F(
                    Var<T,DIM,MEM>& a_prod,
                    const Var<T,DIM,MEM,DIM>& a_left,
                    const Var<T,DIM,MEM>& a_right,
                    T a_scale)
{
  for (int ind0 = 0; ind0 < DIM; ind0++)
    {      
      a_prod(ind0) =0.;
      for (int ind = 0; ind < DIM; ind++)
        {
          a_prod(ind0) +=
            a_left(ind0,ind)*a_right(ind);
        }
      a_prod(ind0) *= a_scale;
    }
}
PROTO_KERNEL_END(f_matVecProd_F, f_matVecProd)

template<typename T, MemType MEM>
void f_matVecProd_(
                   vecType& a_prod,
                   const matType& a_left,
                   const vecType& a_right)
{
  for (int ind0 = 0; ind0 < DIM; ind0++)
    {      
      a_prod[ind0] =0.;
      for (int ind = 0; ind < DIM; ind++)
        {
          a_prod[ind0] +=
            a_left[ind0][ind]*a_right[ind];
        }
    }
}
template<typename T, MemType MEM, unsigned int C, unsigned int D>
PROTO_KERNEL_START
void f_matTransVecProd_F(
                         Var<T,C,MEM>& a_prod,
                         const Var<T,D,MEM,C>& a_left,
                         const Var<T,D,MEM>& a_right,
                         T a_scale)
{
  for (int ind0 = 0; ind0 < C; ind0++)
    {      
      a_prod(ind0) =0.;
      for (int ind = 0; ind < D; ind++)
        {
          a_prod(ind0) +=
            a_left(ind,ind0)*a_right(ind)*a_scale;
        }
    }
}
PROTO_KERNEL_END(f_matTransVecProd_F, f_matTransVecProd)

template<typename T>
void f_matTransVecProd_(
                         vecType& a_prod,
                         const matType& a_left,
                         const vecType& a_right,
                         T a_scale)
{
  for (int ind0 = 0; ind0 < DIM; ind0++)
    {      
      a_prod[ind0] =0.;
      for (int ind = 0; ind < DIM; ind++)
        {
          a_prod[ind0] +=
            a_left[ind][ind0]*a_right[ind]*a_scale;
        }
    }
}

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_matinv3by3_F(
                   Var<T,DIM,MEM,DIM>& a_inverse,
                   const Var<T,DIM,MEM,DIM>& a_matrix)
{
  //cout << "cross-product:" << endl;
  for (int ind0 = 0; ind0 < DIM; ind0++)
    {
      int ind0p1 = (ind0+1)%DIM;
      int ind0p2 = (ind0+2)%DIM;   
      for (int ind1 = 0; ind1 < DIM; ind1++)
        {
          int ind1p1 = (ind1+1)%DIM;
          int ind1p2 = (ind1+2)%DIM;
          a_inverse(ind0,ind1) =
            (a_matrix(ind1p1,ind0p1)*a_matrix(ind1p2,ind0p2) -
             a_matrix(ind1p1,ind0p2)*a_matrix(ind1p2,ind0p1));
        }
    }
  //cout << endl;
  T detMat = 0.;
  for (int ind = 0; ind < DIM; ind++)
    {
      detMat += a_inverse(0,ind)*a_matrix(ind,0);
    }
  for (int ind0 = 0; ind0 < DIM; ind0++)
    {
       for (int ind1 = 0; ind1 < DIM; ind1++)
         {
           a_inverse(ind0,ind1) /= detMat;
         }
    }
}
PROTO_KERNEL_END(f_matinv3by3_F, f_matinv3by3)

template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_Amatrix_F(Point a_pt,
                 Var<T,DIM,MEM,DIM>& a_amatrix,
                 Array<unsigned int,DIM> a_permute,
                 Array<int,DIM> a_sign,
                T a_h,
                T a_offseta,
                T a_offsetb)
{
  T alpha = -.25*M_PI + .5*M_PI*(a_pt[1] + a_offseta)*a_h;
  T tana = tan(alpha);
  T secasq = 1.0 + tana*tana;
           
  T beta =
    -.25*M_PI + .5*M_PI*(a_pt[2] + a_offsetb)*a_h;
  T tanb = tan(beta);
  T secbsq = 1.0 + tanb*tanb;
  T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
  T normsq = norm*norm;
  T norm3 = normsq*norm;
  
  matType amat2;
  amat2[0][0] =  norm;
  amat2[1][0] =  tana*norm;
  amat2[2][0] =  tanb*norm;
  amat2[0][1] = -tana*secasq*norm3*M_PI/2;
  amat2[1][1] =  secasq*secbsq*norm3*M_PI/2;
  amat2[2][1] = -secasq*tana*tanb*norm3*M_PI/2;
  amat2[0][2] = -tanb*secbsq*norm3*M_PI/2;
  amat2[1][2] = -secbsq*tana*tanb*norm3*M_PI/2;
  amat2[2][2] =  secbsq*secasq*norm3*M_PI/2;
  for (int dir0 = 0; dir0 < DIM; dir0++)
    {
      for (int dir1 = 0; dir1 < DIM; dir1++)
        {
          a_amatrix(dir0,dir1) = a_sign[dir0]*amat2[a_permute[dir0]][dir1];
        }
    }        
}

PROTO_KERNEL_END(f_Amatrix_F,f_Amatrix)
template<typename T>
void f_Ainv_(Point a_pt,
             matType& a_amat,
             Array<unsigned int,DIM> a_permute,
             Array<int,DIM> a_sign,
             T a_h,
             T a_offseta,
             T a_offsetb,
             unsigned int a_block)
{
  T alpha = -.25*M_PI + .5*M_PI*(a_pt[1] + a_offseta)*a_h;
  T tana = tan(alpha);
  T secasq = 1.0 + tana*tana;
           
  T beta =
    -.25*M_PI + .5*M_PI*(a_pt[2] + a_offsetb)*a_h;
  T tanb = tan(beta);
  T secbsq = 1.0 + tanb*tanb;
  
  T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
  T normsq = norm*norm;
  T norm3 = normsq*norm;

  matType amat2;
  amat2[0][0] = 1.0/norm;
  amat2[0][1] = tana*norm;
  amat2[0][2] = tanb*norm;
  amat2[1][0] =-tana/(secasq*norm*M_PI/2);
  amat2[1][1] = 1.0/(secasq*norm*M_PI/2);
  amat2[1][2] = 0.0;
  amat2[2][0] = -tanb/(secbsq*norm*M_PI/2);
  amat2[2][1] = 0.0;
  amat2[2][2] = 1.0/(secbsq*norm*M_PI/2);
  for (int dir0 = 0; dir0 < DIM ; dir0++)
    {
      for (int dir1 = 0; dir1 < DIM ; dir1++)
        {
          a_amat[dir0][dir1] = a_sign[dir1]*amat2[dir0][a_permute[dir1]];
        }
    }
}

template<typename T>
void f_detA_A_(Point a_pt,
               matType& a_deta_amat,
               Array<unsigned int,DIM> a_permute,
               Array<int,DIM> a_sign,
               T a_h,
               T a_offseta,
               T a_offsetb)
{
  T alpha = -.25*M_PI + .5*M_PI*(a_pt[1] + a_offseta)*a_h;
  T tana = tan(alpha);
  T secasq = 1.0 + tana*tana;
           
  T beta =
    -.25*M_PI + .5*M_PI*(a_pt[2] + a_offsetb)*a_h;
  T tanb = tan(beta);
  T secbsq = 1.0 + tanb*tanb;
  
  T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
  T normsq = norm*norm;
  T norm3 = normsq*norm;
  T deta = secasq*secbsq*norm3*M_PI*M_PI*.25;
  
   matType amat2;
   amat2[0][0] =  norm;
   amat2[1][0] =  tana*norm;
   amat2[2][0] =  tanb*norm;
   amat2[0][1] = -tana*secasq*norm3*M_PI/2;
   amat2[1][1] =  secasq*secbsq*norm3*M_PI/2;
   amat2[2][1] = -secasq*tana*tanb*norm3*M_PI/2;
   amat2[0][2] = -tanb*secbsq*norm3*M_PI/2;
   amat2[1][2] = -secbsq*tana*tanb*norm3*M_PI/2;
   amat2[2][2] =  secbsq*secasq*norm3*M_PI/2;
   
   for (int dir0 = 0; dir0 < DIM ; dir0++)
     {
       for (int dir1 = 0; dir1 < DIM ; dir1++)
         {
           a_deta_amat[dir0][dir1] = a_sign[dir0]*amat2[a_permute[dir0]][dir1]*deta;
         }
     }
}

template<typename T>
void f_Amatrix_(Point a_pt,
               matType& a_deta_amat,
               Array<unsigned int,DIM> a_permute,
               Array<int,DIM> a_sign,
               T a_h,
               T a_offseta,
               T a_offsetb)
{
  T alpha = -.25*M_PI + .5*M_PI*(a_pt[1] + a_offseta)*a_h;
  T tana = tan(alpha);
  T secasq = 1.0 + tana*tana;
           
  T beta =
    -.25*M_PI + .5*M_PI*(a_pt[2] + a_offsetb)*a_h;
  T tanb = tan(beta);
  T secbsq = 1.0 + tanb*tanb;
  
  T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
  T normsq = norm*norm;
  T norm3 = normsq*norm;
  T deta = secasq*secbsq*norm3;
  
   matType amat2;
   amat2[0][0] =  norm;
   amat2[1][0] =  tana*norm;
   amat2[2][0] =  tanb*norm;
   amat2[0][1] = -tana*secasq*norm3;
   amat2[1][1] =  secasq*secbsq*norm3;
   amat2[2][1] = -secasq*tana*tanb*norm3;
   amat2[0][2] = -tanb*secbsq*norm3;
   amat2[1][2] = -secbsq*tana*tanb*norm3;
   amat2[2][2] =  secbsq*secasq*norm3;
   
   for (int dir0 = 0; dir0 < DIM ; dir0++)
     {
       for (int dir1 = 0; dir1 < DIM ; dir1++)
         {
           a_deta_amat[dir0][dir1] = a_sign[dir0]*amat2[a_permute[dir0]][dir1];
         }
     }
}

template<typename T,MemType MEM>
PROTO_KERNEL_START
void f_detA_F(Point a_pt,
              Var<T,1,MEM>& a_detA,
              T a_h,
              T a_offseta,
              T a_offsetb)
{  
  T alpha = -.25*M_PI + .5*M_PI*(a_pt[1] + a_offseta)*a_h;
  T tana = tan(alpha);
  T secasq = (1.0 + tan(alpha)*tan(alpha));          
  T beta =
    -.25*M_PI + .5*M_PI*(a_pt[2] + a_offsetb)*a_h;
  T tanb = tan(beta);
  T secbsq = (1.0 + tan(beta)*tan(beta));
  T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
  T norm3 = norm*norm*norm;
  a_detA(0) = secasq*secbsq*norm3;
}
PROTO_KERNEL_END(f_detA_F,f_detA)
template<typename T,MemType MEM>
void f_detA_(Point a_pt,
             T& a_detA,
             T a_h,
             T a_offseta,
             T a_offsetb)
{  
  T alpha = -.25*M_PI + .5*M_PI*(a_pt[1] + a_offseta)*a_h;
  T tana = tan(alpha);
  T secasq = (1.0 + tan(alpha)*tan(alpha));          
  T beta =
    -.25*M_PI + .5*M_PI*(a_pt[2] + a_offsetb)*a_h;
  T tanb = tan(beta);
  T secbsq = (1.0 + tan(beta)*tan(beta));
  T norm = 1.0/sqrt(tana*tana + tanb*tanb + 1.0);
  T norm3 = norm*norm*norm;
  a_detA = secasq*secbsq*norm3*(M_PI*M_PI)*.25;
}
#if 0
template<typename T, MemType MEM>
PROTO_KERNEL_START
void f_consToPrim_F(
        Var<T, NUMCOMPS, MEM>&          a_W,
        const Var<T,DIM, MEM,DIM>&      a_Amat,
        const Var<T, NUMCOMPS, MEM>&    a_U,
        double                          a_gamma)
{
    double rho = a_U(0);
    double v2 = 0.0;
    a_W(0) = rho;

    for (int i = 1; i <= DIM; i++)
    {
        double v;
        v = a_U(i) / rho;

        a_W(i) = v;
        v2 += v*v;
    }

    a_W(NUMCOMPS-1) = (a_U(NUMCOMPS-1) - .5 * rho * v2) * (a_gamma - 1.0);
    
}
PROTO_KERNEL_END(f_consToPrim_F, f_consToPrim);
#endif

template<MemType MEM>
void f_cubeSphereMapF(
    Point& a_pt, Var<double,DIM,MEM>& a_X, Var<double,1,MEM>& a_J,
    Array<double, DIM> a_dx, double a_r0, double a_r1,
    int a_block, int a_radialCoord)
{
    // subscript 0 are at corners
    Array<double, DIM> x  = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);

    double r0, xi0, eta0, r, xi, eta;
    switch (a_radialCoord)
    {
        case 0:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[0];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[1];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[2];
            break;
        }
        case 2:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[2];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[0];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[1];
            break;
        }
    }
    
    double X0 = tan(xi0);
    double Y0 = tan(eta0);
    double d0 = sqrt(1+X0*X0+Y0*Y0);
    
    double X = tan(xi);
    double Y = tan(eta);
    double d = sqrt(1+X*X+Y*Y);

    a_J(0) = r0*r0*(1+X0*X0)*(1+Y0*Y0)/(d0*d0*d0);
    switch (a_block)
    {
        case 0:
            a_X(0) = r0*X0/d0;
            a_X(1) = -r0*Y0/d0;
            a_X(2) = -r0/d0;

            break;
        case 1:
            a_X(0) = r0*X0/d0;
            a_X(1) = r0*Y0/d0;
            a_X(2) = r0/d0;
            break;
        case 2:
            a_X(0) = r0*X0/d0;
            a_X(1) = -r0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        case 3:
            a_X(0) = r0/d0;
            a_X(1) = r0*X0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        case 4:
            a_X(0) = -r0*X0/d0;
            a_X(1) = r0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        case 5:
            a_X(0) = -r0/d0;
            a_X(1) = -r0*X0/d0;
            a_X(2) = r0*Y0/d0;
            break;
        default:
            a_X(0) = 0;
            a_X(1) = 0;
            a_X(2) = 0;
            break;
    }
}
PROTO_KERNEL_END(f_cubeSphereMapF, f_cubeSphereMap);

template<MemType MEM>
void f_cubeSphereMapPoint(
    Point& a_pt, vecType& a_X, Array<double, DIM> a_dx, double a_r0, double a_r1,
    int a_block, int a_radialCoord)
{
    // subscript 0 are at corners
    Array<double, DIM> x  = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);

    double r0, xi0, eta0, r, xi, eta;
    switch (a_radialCoord)
    {
        case 0:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[0];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[1];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[2];
            break;
        }
        case 2:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[2];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[0];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[1];
            break;
        }
    }
    
    double X0 = tan(xi0);
    double Y0 = tan(eta0);
    double d0 = sqrt(1+X0*X0+Y0*Y0);
    
    double X = tan(xi);
    double Y = tan(eta);
    double d = sqrt(1+X*X+Y*Y);

    switch (a_block)
    {
        case 0:
            a_X[0] = r0*X0/d0;
            a_X[1] = -r0*Y0/d0;
            a_X[2] = -r0/d0;

            break;
        case 1:
            a_X[0] = r0*X0/d0;
            a_X[1] = r0*Y0/d0;
            a_X[2] = r0/d0;
            break;
        case 2:
            a_X[0] = r0*X0/d0;
            a_X[1] = -r0/d0;
            a_X[2] = r0*Y0/d0;
            break;
        case 3:
            a_X[0] = r0/d0;
            a_X[1] = r0*X0/d0;
            a_X[2] = r0*Y0/d0;
            break;
        case 4:
            a_X[0] = -r0*X0/d0;
            a_X[1] = r0/d0;
            a_X[2] = r0*Y0/d0;
            break;
        case 5:
            a_X[0] = -r0/d0;
            a_X[1] = -r0*X0/d0;
            a_X[2] = r0*Y0/d0;
            break;
        default:
            a_X[0] = 0;
            a_X[1] = 0;
            a_X[2] = 0;
            break;
    }
}

template<MemType MEM>
void f_cubeSphereMap2F(
    Point& a_pt, Var<double,DIM,MEM>& a_X, Var<double,1,MEM>& a_J,
    Var<double,DIM,MEM>& a_n_x_d2n, Var<double,DIM,MEM>& a_d1n_x_n,
    Array<double, DIM> a_dx, double a_r0, double a_r1,
    int a_block, int a_radialCoord)
{
    // subscript 0 are at corners
    Array<double, DIM> x  = a_pt;
    x += 0.5;
    x *= a_dx;
    Array<double, DIM> x0 = x - (0.5*a_dx);

    double r0, xi0, eta0, r, xi, eta;
    switch (a_radialCoord)
    {
        case 0:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[0];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[1];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[2];
            break;
        }
        case 2:
        {
            r0   = a_r0 + (a_r1-a_r0)*x0[2];
            xi0  = -M_PI/4.0 + M_PI/2.0*x0[0];
            eta0 = -M_PI/4.0 + M_PI/2.0*x0[1];
            break;
        }
    }
    
    double X0 = tan(xi0);
    double Y0 = tan(eta0);
    double d0 = sqrt(1+X0*X0+Y0*Y0);
    a_J(0) = r0*r0*(1+X0*X0)*(1+Y0*Y0)/(d0*d0*d0);
    r0 = a_J(0);
    
    Array<double,DIM> cross1; // n_x_d2n
    cross1[0] =  (1 + Y0*Y0)*X0/(d0*d0);    
    cross1[1] = -(1 + Y0*Y0)/(d0*d0);
    cross1[2] = 0.;
    Array<double,DIM> cross2; // d1n_x_n
    cross2[0] =  (1 + X0*X0)*Y0/(d0*d0);
    cross2[1] = 0.;
    cross2[2] = -(1 + X0*X0)/(d0*d0);
    
    double X = tan(xi);
    double Y = tan(eta);
    double d = sqrt(1+X*X+Y*Y);
 
    switch (a_block)
    {
        case 0:
            a_X(0) = r0*X0/d0;
            a_X(1) = -r0*Y0/d0;
            a_X(2) = -r0/d0;
            a_n_x_d2n(0) = cross1[1];
            a_n_x_d2n(1) = -cross1[2];
            a_n_x_d2n(2) = -cross1[0];
            a_d1n_x_n(0) = cross2[1];
            a_d1n_x_n(1) = -cross2[2];
            a_d1n_x_n(2) = -cross2[0];
            break;
        case 1:
            a_X(0) = r0*X0/d0;
            a_X(1) = r0*Y0/d0;
            a_X(2) = r0/d0;
            a_n_x_d2n(0) = cross1[1];
            a_n_x_d2n(1) = cross1[2];
            a_n_x_d2n(2) = cross1[0];
            a_d1n_x_n(0) = cross2[1];
            a_d1n_x_n(1) = cross2[2];
            a_d1n_x_n(2) = cross2[0];
            break;
        case 2:
            a_X(0) = r0*X0/d0;
            a_X(1) = -r0/d0;
            a_X(2) = r0*Y0/d0;
            a_n_x_d2n(0) = cross1[1];
            a_n_x_d2n(1) = -cross1[0];
            a_n_x_d2n(2) = cross1[2];
            a_d1n_x_n(0) = cross2[1];
            a_d1n_x_n(1) = -cross2[0];
            a_d1n_x_n(2) = cross2[2];
            break;
        case 3:
            a_X(0) = r0/d0;
            a_X(1) = r0*X0/d0;
            a_X(2) = r0*Y0/d0;
            a_n_x_d2n(0) = cross1[0];
            a_n_x_d2n(1) = cross1[1];
            a_n_x_d2n(2) = cross1[2];
            a_d1n_x_n(0) = cross2[0];
            a_d1n_x_n(1) = cross2[1];
            a_d1n_x_n(2) = cross2[2];
            break;
        case 4:
            a_X(0) = -r0*X0/d0;
            a_X(1) = r0/d0;
            a_X(2) = r0*Y0/d0;
            a_n_x_d2n(0) = -cross1[1];
            a_n_x_d2n(1) = cross1[0];
            a_n_x_d2n(2) = cross1[2];
            a_d1n_x_n(0) = -cross2[1];
            a_d1n_x_n(1) = cross2[0];
            a_d1n_x_n(2) = cross2[2];
            break;
        case 5:
            a_X(0) = -r0/d0;
            a_X(1) = -r0*X0/d0;
            a_X(2) = r0*Y0/d0;
            a_n_x_d2n(0) = -cross1[0];
            a_n_x_d2n(1) = -cross1[1];
            a_n_x_d2n(2) = cross1[2];
            a_d1n_x_n(0) = -cross2[0];
            a_d1n_x_n(1) = -cross2[1];
            a_d1n_x_n(2) = cross2[2];
            break;
        default:
            a_X(0) = 0;
            a_X(1) = 0;
            a_X(2) = 0;
            break;
    }
}
PROTO_KERNEL_END(f_cubeSphereMap2F, f_cubeSphereMap2);

template<MemType MEM>
PROTO_KERNEL_START
void f_cartesianToSphericalMapF(Point& a_pt, Var<double,DIM,MEM>& a_R,
                                Var<double,1,MEM>& a_J,
        Var<double, DIM, MEM>& a_X, int a_block, int a_radialCoord)
{
    double x0 = a_X(0);
    double y0 = a_X(1);
    double z0 = a_X(2);

    double x, y, z;
    switch (a_block)
    {
        case 0:
        {
            x = -z0;
            y = y0;
            z = x0;
            break;
        }
        case 1:
        {
            x = z0;
            y = y0;
            z = -x0;
            break;
        }
        case 2:
        {
            x = -y0;
            y = x0;
            z = z0;
            break;
        }
        case 3:
        {
            x = x0;
            y = y0;
            z = z0;
            break;
        }
        case 4:
        {
            x = y0;
            y = -x0;
            z = z0;
            break;
        }
        case 5:
        {
            x = -x0;
            y = -y0;
            z = z0;
            break;
        }
    }

    double R = sqrt(x*x + y*y + z*z);
    double r = sqrt(x*x + y*y);
    double phi = atan2(r,z);
    double theta = atan2(y,x);
    // note that this J is at corners
    a_J(0) = r*r*sin(phi);
    a_R(a_radialCoord) = R;
    a_R((a_radialCoord + 1) % 3) = theta;
    a_R((a_radialCoord + 2) % 3) = phi;
}
PROTO_KERNEL_END(f_cartesianToSphericalMapF, f_cartesianToSphericalMap);
#endif // end include guard.
