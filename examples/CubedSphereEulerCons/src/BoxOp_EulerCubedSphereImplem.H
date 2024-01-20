template<typename T,MemType MEM>
void
primToFlux(
           State& a_physFlux,
           const State& a_WFace,
           const State& a_WFaceBar,
           const BoxData<T,1,MEM>& a_adjDdir,
           T a_h,
           T a_gamma,
           uint a_block,
           uint a_centering)
{
  // Compute fluxes at face centers.
  PR_TIMERS("primToFlux");
  Stencil<T> deconvolveStencilRadial = (1.0*Shift(Point::Zeros())
                                    - (1.0/24.0)*Stencil<T>::Derivative(2,0,2));
  // Deconvolve W to obtain point values.
  BoxData<T,NUMCOMPS,MEM> WFacePoint = _deconvolveFace(a_WFace,a_WFaceBar,a_centering);  
  BoxData<T,1,MEM> kePoint = kineticEnergyW<T,MEM>(WFacePoint,a_h,a_block,a_centering);
 
  BoxData<T,1,MEM> adjDdirPoint;
  if (a_centering != 0)
    {
      adjDdirPoint = deconvolveStencilRadial(a_adjDdir);
    }
  else
    {
      adjDdirPoint.define(a_adjDdir.box());
      a_adjDdir.copyTo(adjDdirPoint);
    }
  State fluxPoint = forall<T,NUMCOMPS>
    ([ ] PROTO_LAMBDA
     (
      Var<T,NUMCOMPS>& a_flux,
      const Var<T,NUMCOMPS>& a_wface,      
      const Var<T,1>& a_ke,
      const Var<T>& a_adjddir,
      const T a_gamma,
      uint a_block,
      uint a_centering)
     {
       T velnorm = a_wface(a_centering+1);
       a_flux(0) = a_wface(0)*velnorm*a_adjddir(0);
       T veltan1 = a_wface((a_centering+1)%DIM + 1);
       T veltan2 = a_wface((a_centering+2)%DIM + 1);
       a_flux(a_centering+1) = a_flux(0)*velnorm;
       a_flux((a_centering+1)%DIM + 1) = a_flux(0)*veltan1;
       a_flux((a_centering+2)%DIM + 1) = a_flux(0)*veltan2;
       a_flux(NUMCOMPS - 1) =
       (a_ke(0) + a_wface(NUMCOMPS-1)*a_gamma/(a_gamma - 1.0))*a_adjddir(0)*velnorm;
     },WFacePoint,kePoint,adjDdirPoint,a_gamma,a_block,a_centering);
   if ((a_block == BLOCKTEST) && a_centering == 0)
    {
      HDF5Handler h5;
      h5.writePatch(a_h,kePoint,"kePoint" + to_string(a_centering));
      h5.writePatch(a_h,WFacePoint,"WPointPF" + to_string(a_centering));          
    }
  // Compute Average fluxes from pointwise values.
  a_physFlux = (1.0*Shift(Point::Zeros()) +
                (1.0/24.0)*Stencil<T>::LaplacianFace(a_centering))(fluxPoint);
}

