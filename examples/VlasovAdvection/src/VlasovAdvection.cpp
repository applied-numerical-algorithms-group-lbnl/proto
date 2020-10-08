#include "VlasovAdvection.H"
#include "Proto.H"
#include "Proto_BoxData.H"

PROTO_KERNEL_START
void computeNormalVelocityF(Var<double,1>& a_veln,
                            Var<double,DIM>& a_metric,
                            Var<double,DIM>& a_vel)
{
    a_veln(0)=0.0;
    for(int k=0; k<DIM; k++)
        a_veln(0)+=a_metric(k)*a_vel(k);
}
PROTO_KERNEL_END(computeNormalVelocityF,computeNormalVelocity)


PROTO_KERNEL_START
void computePhiTimesVelF(Var<double,DIM>& a_phivel,
                         Var<double,1>& a_phi,
                         Var<double,DIM>& a_vel)
{
    for(int k=0; k<DIM; k++)
        a_phivel(k)=a_phi(0)*a_vel(k);
}
PROTO_KERNEL_END(computePhiTimesVelF,computePhiTimesVel)

PROTO_KERNEL_START
void add4thOrderCorrectionVecF(Var<double,DIM>& a_data,
                               Var<double,DIM>& a_correction,
                               const double& a_dx)
{
    double scale=a_dx*a_dx/12.0;
    for(int k=0; k<DIM; k++)
        a_data(k)+=scale*a_correction(k);
}
PROTO_KERNEL_END(add4thOrderCorrectionVecF,add4thOrderCorrectionVec)

PROTO_KERNEL_START
void add4thOrderCorrectionScalarF(Var<double,1>& a_data,
                                  Var<double,1>& a_correction,
                                  const double& a_dx)
{
    double scale=a_dx*a_dx/12.0;
    a_data(0)+=scale*a_correction(0);
}
PROTO_KERNEL_END(add4thOrderCorrectionScalarF,add4thOrderCorrectionScalar)

void setCellCenteredVelocity(BoxData<double,DIM>& vel,
                             Box& patch,
                             const int idir)
{
    //TODO: Fill in
}

void setFaceAveragedMetrics(BoxData<double,DIM>& metrics,
                            Box& patch,
                            const int idir)
{
    //TODO: Fill in
}

void BWENO(BoxData<double,1>& phi_faceave,
           BoxData<double,1>& normal_vel,
           BoxData<double,1>& phi_cellave)
{
    //TODO: Fill in
}



/////////////////////////////////////////////////////////////////////////
//TODO: Velocities and metrics will be generated patch-by-patch
//Transcribe velocities, metric computations from Fortran to Proto
//Document algorithm
//Move loop over idir into data iterator loop
void computeRHS(LevelBoxData<double,1>& a_rhs,
                const LevelBoxData<double,1>& a_phi,
                const double a_dx)
{
    Stencil<double> grad=Stencil<double>::Derivative(1,0,2);
    for(int idir=1; idir<DIM; idir++)
    {
        grad+=Stencil<double>::Derivative(1,idir,2);
    }

    double div_flux_scale=-1.0/a_dx;
    a_rhs.define(a_phi.getDBL(),Point::Zero());
    a_rhs.setToZero();
    a_phi.exchange();
    for(DataIterator dit=a_phi.begin(); *dit!=dit.end(); ++dit)
    {
        Box& patch_box=a_phi[*dit].box();
        BoxData<double,1>& rhs_bd=rhs[*dit];
        for(int idir=0; idir<DIM; idir++)
        {
            Stencil<double> lap_transv=LaplacianFace(idir);
            //Gradient transverse to idir
            Stencil<double> grad_transv=grad-Stencil<double>::Derivative(1,idir,2);

            //////////////////////////////////////////////////////////////////////////////////
            //Setup the velocity and metrics
            BoxData<double,DIM> velfc_bd;
            setCellCenteredVelocity(velfc_bd,patch_box,idir);
            BoxData<double,DIM> metrics;
            setFaceAveragedMetrics(metrics,patch_box,idir);

            //////////////////////////////////////////////////////////////////////////////////
            //Compute face averages of velocity
            //Box that's shrunk by 1 in the transverse direction.
            Box fa_box=velfc_bd.box().grow(-1).grow(idir,1);
            BoxData<double,DIM> velfa_bd(fa_box);
            velfa_bd|=velfc_bd;
            velfa_bd+=lap_transv(velfc_bd,a_dx*a_dx/24.0);

            //////////////////////////////////////////////////////////////////////////////////
            //Compute face averages of phi
            BoxData<double,1> normal_vel(velfa_bd.box());
            forallInPlace(normal_vel,metrics,velfa_bd);
            //Can I use the dit from a_vel for a_phi when it has a totally different centering?
            //Is this the correct way to set up the box?
            //Actually need to reduce this box by 2 to account for ghost cells in a_phi...
            BoxData<double,1> phifa_bd(patch_box.growHi(idir,1));
            BWENO(phifa_bd,normal_vel,a_phi[*dit]);

            ////////////////////////////////////////////////////////////////////////////////////
            //Compute 4th order face average of phi*vel
            BoxData<double,DIM> pvfa_bd(phifa_bd.box());
            forallInPlace(computePhiTimesVel,pvfa_bd,phifa_bd,velfa_bd);
            BoxData<double,DIM> gradpvfa_bd=grad_transv(pvfa_bd); //Stored for next step
            BoxData<double,1> gradphifa_bd=grad_transv(phifa_bd);
            BoxData<double,DIM> gradvelfa_bd=grad_transv(velfa_bd);
            BoxData<double,DIM> pvcorrection_bd(gradphifa_bd.box());
            forallInPlace(computePhiTimesVel,pvcorrection_bd,gradphifa_bd,gradvelfa_bd);
            forallInPlace(add4thOrderCorrectionVec,pvfa_bd,pvcorrection_bd,a_dx);

            ////////////////////////////////////////////////////////////////////////////////////
            //Compute final flux
            BoxData<double,DIM> flux(pvfa_bd.box());
            BoxData<double,DIM> gradMetric_bd=grad_transv(metrics);
            BoxData<double,1> flux_correction_bd(gradMetric_bd.box());
            forallInPlace(computeNormalVelocity,flux_correction_bd,gradMetric_bd,gradpvfa_bd);
            forallInPlace(computeNormalVelocity,flux,metrics,pvfa_bd);
            forallInPlace(add4thOrderCorrectionScalar,flux,flux_correction_bd,a_dx);
            rhs_bd+=Stencil<double>::FluxDivergence(idir)(flux,div_flux_scale);
        }
    }
}
