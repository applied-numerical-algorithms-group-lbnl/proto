#ifndef _VLASOV_ADVECTION_H
#define _VLASOV_ADVECTION_H

class VlasovAdvection
{
public:
    VlasovAdvection(vector<vector<BoxData<double>>>& N_mat_);
    /**
     * Computes -div(N^Tu\phi)_{X_i}, where N^T is the mapping from configuration to physical space,
     * u is the velocity, and \phi is the distribution.
     */
    void computeRHS(BoxData<double>& rhs, BoxData<double>& phi_cell, vector<BoxData<double>>& vel_cent);
    void computeVelFaceAve(vector<BoxData<double>>& vel_ave, vector<BoxData<double>>& vel_cent, int idir);
    void computePhiFaceAve(BoxData<double>& phi_face, BoxData<double>& phi_cell, vector<BoxData<double>>& vel_ave, int idir);
    void computeVelPhiFaceAve(vector<BoxData<double>>& phivel, BoxData<double>& phi, vector<BoxData<double>>& vel, int idir);
    void computeProductFaceAve(vector<BoxData<double>>& flux, vector<BoxData<double>>& phivel, int idir);
    void computeTransverseGradient(BoxData<double>& data, int idir);

private:
//Mapping from configuration to physical space
    vector<vector<BoxData<double>>> N_mat;
    vector<vector<BoxData<double>>> trans_grad_uphi;

    void ComputeWENOSmoothnessFactors(/*fill in*/);
    void ComputeWENOFlux(/*fill in*/);
};

#endif
