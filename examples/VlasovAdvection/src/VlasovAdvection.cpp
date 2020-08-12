#include "VlasovAdvection.H"

VlasovAdvection::VlasovAdvection(vector<vector<BoxData<double>>>& N_mat_,
                                 const float& dx_) {
    N_mat=N_mat_;
    dx=dx_;
}

void VlasovAdvection::computeRHS(BoxData<double>& rhs, BoxData<double>& phi_cell, vector<BoxData<double>>& vel_cent) {
    //Set up rhs domain, values, initialize to 0.0
//Compute flux in each dir and accumulate difference into rhs
    for(int idir=0; idir<vel_cent.size(); idir++) {
        vector<BoxData<double>> vel=vel_cent; //May not be able to assign BoxData like this
        for(int k=0; k<vel_cent.size(); k++) {
            //Compute the Laplacian
            //Add to the vel (will this keep the same dimensions as vel_cent and not affect the ghost faces)?
        }
        BoxData<double> Nvel(Nvel[idir][0].box());
        for(int k=0; k<vel.size(); k++) {
            Nvel+=Nvel_2nd[idir][k]*vel[k]; //Probably not the right way to multiply-accumulate
        }
        BoxData<double> phi_face(Nvel.box());
        WENOFaceAverageFromFaceCell(phi_face,Nvel,phi_cell,idir);
        vector<BoxData<double>> phi_vel_face(vel.size());
        for(int k=0; k<vel.size(); k++) {
            phi_vel_face[k]=phi_face*vel[k];
        }
        vector<BoxData<double>> grad_phivel(vel.size());
        for(int k=0; k<vel.size(); k++) {
            grad_phivel[k]=TransverseGradient(phi_vel_face[k],idir);
            BoxData<double> grad_phi=TransverseGradient(phi_face,idir);
            BoxData<double> grad_velk=TransverseGradient(vel[k],idir);
            BoxData<double> correction=grad_phi*grad_velk;
            correction*=(dx*dx/12);
            phi_vel_face+=correction;
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



