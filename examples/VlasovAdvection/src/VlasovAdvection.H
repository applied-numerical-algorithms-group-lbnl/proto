#ifndef _VLASOV_ADVECTION_H
#define _VLASOV_ADVECTION_H

class VlasovAdvection
{
public:
    VlasovAdvection(vector<vector<BoxData<double>>>& N_mat);
    void computeRHS(BoxData<double>& rhs, BoxData<double>& phi_cell, vector<BoxData<double>>& vel_cent);
    void computeVelFaceAve(BoxData<double>& vel_ave, BoxData<double>& vel_cent, int idir);
    void computePhiFaceAve(BoxData<double>& phi_face, BoxData<double>& phi_cell);
    void computeProductFaceAve(BoxData<double>& phivel, BoxData<double>& phi, BoxData<double>& vel, int idir);
    void computeProductFaceAve(vector<BoxData<double>>& flux, vector<BoxData<double>>& phivel);
    void computeTransverseGradient(BoxData<double>& data, int idir);

private:
//Mapping from configuration to physical space
    vector<vector<BoxData<double>>>& N_mat;

    void ComputeWENOSmoothnessFactors(/*fill in*/);
    void ComputeWENOFlux(/*fill in*/);
};

#endif
