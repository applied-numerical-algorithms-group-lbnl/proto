#include "VlasovAdvection.H"

VlasovAdvection::VlasovAdvection(vector<vector<BoxData<double>>>& N_mat_,
                                 const float& dx_) {
    N_mat=N_mat_;
    dx=dx_;
}

BoxData<double> VlasovAdvection::FaceCenterToFaceAverage(BoxData<double>& face_center,
                                                         const int idir)
{
    //May want to make this static member depending on the cost for regenerating stencils
    Stencil<double> Lap2nd=Stencil<double>::Laplacian();
    //TODO: Put code here. Box needs to be shrunk in all directions except idir.
    //See the corresponding function in HOAMRTools
    BoxData<double> face_average(face_center.box());
    return face_average;
}

void VlasovAdvection::WENOFaceAverageFromFaceCell(BoxData<double>& phi_face,
                                                  BoxData<double>& vel,
                                                  BoxData<double>& phi_cell,
                                                  const int idir)
{
    //Move over code for computing WENO flux
}

BoxData<double> VlasovAdvection::TransverseGradient(BoxData<double>& data,
                                                    const int idir)
{
    BoxData<double> grad;
    //Add code for computing the transverse gradient from stencils
}

void VlasovAdvection::computeRHS(BoxData<double>& rhs,
                                 BoxData<double>& phi_cell,
                                 vector<BoxData<double>>& vel_cent,
                                 int phi_ghost,
                                 int vel_ghost)
{
    rhs.define(phi_cell.box().grow(-phi_ghost));
    rhs.setVal(0.0);
//Compute flux in each dir and accumulate difference into rhs
    for(int idir=0; idir<vel_cent.size(); idir++) {
        vector<BoxData<double>> vel(vel_cent.size()); //May not be able to assign BoxData like this
        for(int k=0; k<vel_cent.size(); k++)
            vel[k]=FaceCenterToFaceAverage(vel_cent,idir);

        BoxData<double> Nvel(Nvel[idir][0].box());
        Nvel.setVal(0.0);
        for(int k=0; k<vel.size(); k++) {
            Nvel+=Nvel_2nd[idir][k]*vel[k]; //Probably not the right way to multiply-accumulate
        }
        BoxData<double> phi_face(Nvel.box());
        WENOFaceAverageFromFaceCell(phi_face,Nvel,phi_cell,idir);

        vector<BoxData<double>> phi_vel_face(vel.size());
        vector<BoxData<double>> grad_phivel(vel.size()); //Compute and store for computing flux below
        for(int k=0; k<vel.size(); k++) {
            phi_vel_face[k]=phi_face*vel[k];
            grad_phivel[k]=TransverseGradient(phi_vel_face[k],idir);
            BoxData<double> grad_phi=TransverseGradient(phi_face,idir);
            BoxData<double> grad_velk=TransverseGradient(vel[k],idir);
            BoxData<double> correction=grad_phi*grad_velk;
            correction*=(dx*dx/12);
            phi_vel_face[k]+=correction;
        }

        BoxData<double> flux(vel[0].box());
        flux.setVal(0.0);
        for(int k=0; k<vel_cent.size(); k++) {
            BoxData<double> Nphivel(phi_vel_face[k].box());
            Nphivel=Nvel[idir][k]*phi_vel_face[k];
            BoxData<double> grad_N=TransverseGradient(Nvel[idir][k],idir);
            BoxData<double> correction=grad_N*grad_phivel[k];
            correction*=(dx*dx/12.0);
            Nphivel+=correction;
            flux+=Nphivel;
        }
        
        //TODO: stencil to compute difference of fluxes in idir and put into rhs
    }
    
}



