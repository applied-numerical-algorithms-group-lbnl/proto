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

private:
//Mapping from configuration to physical space
    vector<vector<BoxData<double>>> N_mat;
    BoxData<double> FaceCenterToFaceAverage(BoxData<double>& face_center,const int idir);
    BoxData<double> TransverseGradient(BoxData<double>& data, const int idir);
    void WENOFaceAverageFromFaceCell(BoxData<double>& phi_face,
                                     BoxData<double>& vel,
                                     BoxData<double>& phi_cell,
                                     const int idir);
};

#endif
