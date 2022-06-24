
template<class T, unsigned int CL,unsigned char DL, unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E> faceMatrixProductAB(
                     const BoxData<T,CL,MEM,DL,E>& a_ql4,
                     const BoxData<T,CR,MEM,DR,E>& a_qr4,
                     const BoxData<T,CL,MEM,DL,E>& a_ql2,
                     const BoxData<T,CR,MEM,DR,E>& a_qr2,
                     int a_dir)
{
  Box boxout=a_ql4.box()*a_qr4.box()*a_ql2.box()*a_ql2.box();
  BoxData<T,CL,MEM,DR,E> qprod(boxout);
  qprod.setToZero();
  for (unsigned int i = 0;i < CL,i++)
    {
      for (unsigned int j = 0;j < DR,j++)
        {
          for (unsigned int k = 0; k < CR;k++)
            {
              auto qprodComp = slice(qprod,i,j);
              auto ql4Comp = slice(a_ql4,i,k);
              auto qr4Comp = slice(a_qr4,k,j);
              auto ql2Comp = slice(a_ql2,i,k);
              auto qr2Comp = slice(a_qr2,k,j);
              qprodComp+=faceProduct(ql4Comp,qr4Comp,ql2Comp,qr2Comp,a_dir);
            }           
        }
    }
  return qprod;
}
template<class T, unsigned int CL, unsigned int CR, unsigned char DL, unsigned char DR, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E> faceMatrixProductATB(
                     const BoxData<T,CL,MEM,DL,E>& a_ql4,
                     const BoxData<T,CR,MEM,DR,E>& a_qr4,
                     const BoxData<T,CL,MEM,DL,E>& a_ql2,
                     const BoxData<T,CR,MEM,DR,E>& a_qr2,
                     int a_dir)
{
  Box boxout=a_ql4.box()*a_qr4.box()*a_ql2.box()*a_ql2.box();
  BoxData<T,DL,MEM,DR,E> qprod(boxout);
  qprod.setToZero();
  for (unsigned int i = 0;i < DL,i++)
    {
      for (unsigned int j = 0;j < DR,j++)
        {
          for (unsigned int k = 0; k < CR;k++)
            {
              auto qprodComp = slice(qprod,i,j);
              auto ql4Comp = slice(a_ql4,k,i);
              auto qr4Comp = slice(a_qr4,k,j);
              auto ql2Comp = slice(a_ql2,k,i);
              auto qr2Comp = slice(a_qr2,k,j);
              qprodComp+=faceProduct(ql4Comp,qr4Comp,ql2Comp,qr2Comp,a_dir);
            }           
        }
    }
  return qprod;
}
template<class T, unsigned int C,unsigned char D, MemType MEM>
BoxData<T,C,MEM,D,E> faceMatrixProductFlux(
                     const BoxData<T,C,MEM,D>& a_ql4,
                     const BoxData<T,D,MEM>& a_qr4,
                     const BoxData<T,C,MEM,D>& a_ql2,
                     const BoxData<T,D,MEM>& a_qr2,
                           int a_dir)
{
  Box boxout=a_ql4.box()*a_qr4.box()*a_ql2.box()*a_ql2.box();
  BoxData<T,1,MEM> qprod(boxout);
  qprod.setToZero();
    for (unsigned int k = 0; k < D;k++)
            {
              auto ql4Comp = slice(a_ql4,k);
              auto qr4Comp = slice(a_qr4,a_dir,k);
              auto ql2Comp = slice(a_ql2,k);
              auto qr2Comp = slice(a_qr2,);
              qprod+=faceProduct(ql4Comp,qr4Comp,ql2Comp,qr2Comp,a_dir);
            }
        }
    }
  return qprod;
}
// 4th-order approximation of (1/ql)*qr where ql is a scalar and qr is matrix-valued.
template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E> faceMatrixQuotient(
                     const BoxData<T,1,MEM,1.1>& a_ql4,
                     const BoxData<T,C,MEM,D,E>& a_qr4,
                     const BoxData<T,1,MEM,1,1>& a_ql2,
                     const BoxData<T,C,MEM,D,E>& a_qr2,
                           int a_dir)
{
  Box boxout=a_ql4.box()*a_qr4.box()*a_ql2.box()*a_ql2.box();
  BoxData<T,C,MEM,D,E> qq(boxout);
  qq.setToZero();
  for (unsigned int i = 0;i < C,i++)
    {
      for (unsigned int j = 0;j < D,j++)
        {
          for (unsigned int k = 0;k < E,k++)
            {
              auto qqComp = slice(qq,i,j,k);
              auto ql4Comp = slice(a_qr4,i,j,k);
              auto ql2Comp = slice(a_qr2,i,j,k);
              qqComp += faceQuotient(a_ql4,a_qr4Comp,a_ql2,a_qr2Comp,a_dir);
            }
        }
    }
  return qq;
}
template<class T, unsigned int C,unsigned char DL, unsigned int CR, unsigned char D, unsigned char E, MemType MEM>
BoxData<T,C,MEM,D,E> matrixProduct2(,
                     const BoxData<T,C,MEM,DL,E>& a_ql,
                     const BoxData<T,CR,MEM,D,E>& a_qr,
                           int a_dir)
{
  Box boxout=a_ql.box()*a_qr.box();
  BoxData<T,C,MEM,D,E> qprod(boxout);
  qprod.setToZero();
  for (unsigned int i = 0;i < C,i++)
    {
      for (unsigned int j = 0;j < D,j++)
        {
          for (unsigned int k = 0; k < CR;k++)
            {
              auto qprodComp = slice(qprod,i,j);
              auto qlComp = slice(a_ql,i,k);
              auto qrComp = slice(a_qr,k,j);
              forallInPlace([](Var<T>& a_q,Var<T>& a_q1,Var<T>& a_q2)
                            {qprod(0) += ql(0)*qr(0);},qprodComp,qlComp,qrComp);
            }
        }
    }
  return qprod;
}
template <class T,Memtype MEM>
BoxData<T,DIM,MEM> cofactor(
                            const BoxData<T,DIM,MEM>& a_X,
                            int a_dir)
{
  PR_assert((dir > -1) && (dir < DIM));
  
#ifeq DIM 2
  dir2 = (a_dir+1)%2;
  T sign = 1;
  if (dir == 1) sign = -1;
  BoxData<T,DIM,MEM> NT = (sign*Stencil<T>::divergence(dir))(a_X);
  return NT;
#endif
#ifeq DIM 3
  // NT = d/(dxi_d')(X \cross dx/dxi_d'')/2 -  d/(dxi_d'')(X \cross dx/dxi_d')/2
  //    = d/dxi_d'(X) \cross d/dxi_d''(X)
  // d' = a_dir+1 mod DIM, d''= a_dir+2 mod DIM

  array<BoxData<T,1,MEM>,DIM> NTarray;
  BoxData<T,DIM,MEM> NT;
  dp = (a_dir + 1)%DIM;
  dp2 = (a_dir + 2)%DIM;
  BoxData<T,DIM,MEM> dXdp = Stencil<T>::divergence(dp)(X);
  BoxData<T,DIM,MEM> dXdp2 = Stencil<T>::divergence(dp2)(X);
  // Loop over components of NT to compute cross-products.
  for (int d = 0; d < DIM;d++)
    {
      auto NTcomp = slice(NT,d);
      auto Xcomp = slice(X,(d+1)%DIM);
      auto dXcompp = slice(dXdp,(d+2)%DIM);
      auto dXcompp2 = slice(dXdp2,(d+2)%DIM);
      auto NT1 = Stencil<T>::divergence(dp1)
        (Product1D(Xcomp,dxCompp2,Xcomp,dXcompp2,dp2));
      auto NT2 = Stencil<T>::divergence(dp2)
        (Product1D(Xcomp,dxCompp,Xcomp,dXcompp,dp));
      // We compute the cofactor in NTarray, then store in NT, to use
      // box inference to get the correct domain.
      auto NTarray[d] = forall([](Var<T>& a_q1,Var<T>& a_q2)
                               {Var<T>& retval(0) = (a_q1(0) - a_q2(0))/2;
                                return retval;},NT1,NT2);
    }
  NT = forall([](Var<T>& a_q0,Var<T>& a_q1,Var<T>& a_q2)
              {Var<T,DIM>& retval;
               retval(0) = a_q0(0);
               retval(1) = a_q1(0);
               retval(2) = a_q2(0);
               return retval;
              };NTarray[0],NTarray[1],NTarray[2]);
  return NT;
#endif
}
}
template <class T,MemType MEM>
BoxData<T,MEM> jacobian(
                    const BoxData<T,DIM,MEM>& a_X,
                    const array<BoxData<T,DIM,MEM>& a_NT
                    )
{
  
}

template<typename T>
Stencil<T>
Stencil<T>::FaceToCell(int a_dir,int a_order)
{
  PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
            "Stencil::FaceToCell(int a_dir, int a_order) invalid for a_dir = %i.\
            a_dir should be in [0,DIM=%i).",a_dir, DIM);
  Stencil<T> ret;
  switch (a_order)
    {
    case 2:
      {
        T coefs[2] = {1,1};
        ret =  coefs[1]*Shift::Basis(a_dir,1) + coefs[0]*Shift::Zeros();
        ret *= (T)(1/2.0);
        break;
      }
    case 4:
      {
        T coefs[4] = {-1,13,13,-1};
        ret = coefs[0]*Shift::Basis(a_dir,-1) + coefs[1]*Shift::Shift::Zeros()
          +coefs[2]*Shift::Basis(a_dir,1) + coefs[3]*Shift::Basis(a_dir,2);
        ret *= (T)(1/24.0);
        break;
      }
    default:
      {
        std::abort();
        break;
      }
    }
  return ret;
}
template<typename T>
Stencil<T>
Stencil<T>::CornerToFace(int a_dir,int a_order)
{
  PROTO_ASSERT((a_dir >= 0) && (a_dir < DIM),
            "Stencil::FaceToCell(int a_dir, int a_order) invalid for a_dir = %i.\
            a_dir should be in [0,DIM=%i).",a_dir, DIM);
  Stencil<T> ret;
  Stencil<T> sten2=(1.0)*Shift::Zeros();
  for (d = a_dir+1 ; d < a_dir + DIM;d++)
    {
      dmod = d%DIM;
      sten2 = sten2*(.5*Shift::Zeros() + .5*Shift::Basis(dmod,1));
    }
  switch (a_order)
    {
    case 2:
      {
        ret = sten2;
        break;
      }
    case 4:
      {
        ret = sten2 - (1./12.0)*sten2*Stencil<T>::LaplacianFace(a_dir,2)               
      }
    default:
      {
        std::abort();
        break;
      }
    }
  return ret;
}
template<typename T>
Stencil<T>
Stencil<T>::CornerToCell(int a_order)
{
  Stencil<T> ret;
  Stencil<T> sten2=(1.0)*Shift::Zeros();
  for (int dir = 0; dir < DIM; dir++)
    {
      sten2 = sten2*(.5*Shift::Zeros() + .5*Shift::Basis(dir,1));
    }
  switch (a_order)
    {
    case 2:
      {
        ret = sten2;
        break;
      }
    case 4:
      {
        ret = sten2 - (1.0/12.0)*sten2*Stencil<T>::Laplacian();
        break;
      }
    default:
      {
        std::abort();
        break;
      }
    }
  return ret;
}
