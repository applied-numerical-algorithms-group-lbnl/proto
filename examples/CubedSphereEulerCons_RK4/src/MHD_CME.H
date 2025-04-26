#include "MHD_Constants.H"
#include <iomanip>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <time.h>

template <typename T, MemType MEM>
PROTO_KERNEL_START
void define_CMEF(const Point& a_pt,
                Var<T, NUMCOMPS, MEM> &a_CME,
                 Var<T, DIM, MEM> &a_x)
{   

    double x = a_x(0);
    double y = a_x(1);
    double z = a_x(2);

    double r_dom, theta_dom, phi_dom;

    r_dom = sqrt(x*x + y*y + z*z);
    theta_dom = acos(z/r_dom);
    phi_dom = atan2(y,x);

    double R_t = ParseInputs::get_CME_r1(); // In R_Sun
    double R_p = ParseInputs::get_CME_r0(); // In R_Sun
    double theta_HW = ParseInputs::get_CME_halfangle(); // In degrees
    double T0  = 1.0*ParseInputs::get_CME_pol_flux_control(); // Twist/AU
    double n = ParseInputs::get_CME_FRiED_n() ;  // Flattening coefficient
    double d_theta = 0.000001;  // Check if sufficiently small
    double shift = ParseInputs::get_CME_shift(); // In R_Sun
    double epsilon1 = 0.000001;
        
        
    double lat = -ParseInputs::get_CME_lat()*c_PI/180.0;
    double lon = ParseInputs::get_CME_lon()*c_PI/180.0;
    double tilt = ParseInputs::get_CME_tilt()*c_PI/180.0;

    

    double rot_x[9], rot_y[9], rot_z[9], rot_x_inv[9], rot_y_inv[9], rot_z_inv[9];
    
    rot_y[0] = cos(lat);
    rot_y[1] = 0.;
    rot_y[2] = -sin(lat);
    rot_y[3] = 0.;
    rot_y[4] = 1.;
    rot_y[5] = 0.;
    rot_y[6] = sin(lat);
    rot_y[7] = 0.;
    rot_y[8] = cos(lat);

    rot_x[0] = 1.;
    rot_x[1] = 0.;
    rot_x[2] = 0.;
    rot_x[3] = 0.;
    rot_x[4] = cos(tilt);
    rot_x[5] = sin(tilt);
    rot_x[6] = 0.;
    rot_x[7] = -sin(tilt);
    rot_x[8] = cos(tilt);

    rot_z[0] = cos(lon);
    rot_z[1] = sin(lon);
    rot_z[2] = 0.;
    rot_z[3] = -sin(lon);
    rot_z[4] = cos(lon);
    rot_z[5] = 0.;
    rot_z[6] = 0.;
    rot_z[7] = 0.;
    rot_z[8] = 1.;

    rot_y_inv[0] = cos(lat);
    rot_y_inv[1] = 0.;
    rot_y_inv[2] = sin(lat);
    rot_y_inv[3] = 0.;
    rot_y_inv[4] = 1.;
    rot_y_inv[5] = 0.;
    rot_y_inv[6] = -sin(lat);
    rot_y_inv[7] = 0.;
    rot_y_inv[8] = cos(lat);

    rot_x_inv[0] = 1.;
    rot_x_inv[1] = 0.;
    rot_x_inv[2] = 0.;
    rot_x_inv[3] = 0.;
    rot_x_inv[4] = cos(tilt);
    rot_x_inv[5] = -sin(tilt);
    rot_x_inv[6] = 0.;
    rot_x_inv[7] = sin(tilt);
    rot_x_inv[8] = cos(tilt);

    rot_z_inv[0] = cos(lon);
    rot_z_inv[1] = -sin(lon);
    rot_z_inv[2] = 0.;
    rot_z_inv[3] = sin(lon);
    rot_z_inv[4] = cos(lon);
    rot_z_inv[5] = 0.;
    rot_z_inv[6] = 0.;
    rot_z_inv[7] = 0.;
    rot_z_inv[8] = 1.;
    
    R_t = R_t*c_SR; // In cm
    R_p = R_p*c_SR; // In cm
    double R0_Vandas = ParseInputs::get_CME_Vandas_R0()*c_SR; // In cm
    
    shift = shift*c_SR; // In cm
    T0 = T0/1.496e13; // Twist/cm
    double theta_HW_rad = theta_HW*c_PI/180.;
    double d_theta_rad = d_theta*c_PI/180.;
    
    //B0 = T0*T0/log(T0*T0*R_p*R_p + 1.0)
    double B0 = ParseInputs::get_CME_a1();

    double x_SR, y_SR, z_SR, x_cm, y_cm, z_cm, rx, ry;
    double theta, theta_1, theta_2, error, theta_center, f_theta_1, f_theta_2, f_theta_center, theta_rad;
    double D, r1, r2, R, n_vec_i, n_vec_j, n_vec_mag, n_hat_i, n_hat_j, r_dash, x_dash, y_dash, z_dash, phi;
    double B0_local, B_x_dash, B_y_dash, B_z_dash, t, B_x, B_y, B_z, theta_of_xyz;
    double x_temp, y_temp, z_temp, x_temp2, y_temp2, z_temp2;
    double x_SR_dom, y_SR_dom, z_SR_dom, x_cm_dom, y_cm_dom, z_cm_dom, B_x_dom, B_y_dom, B_z_dom;
    double arg, temp1, temp2;
    double r_1, r_2, R_big, sqBTotal, e0TD;
    double a_hat_i, a_hat_j, a_hat_k, b_hat_i, b_hat_j, b_hat_k, D1, V_rad, V_exp, V_x, V_y, V_z, V_x_dom, V_y_dom, V_z_dom;

    double half = 0.5;
    double d_1_4PI = 1.0/(4.0*c_PI);

    x_cm_dom=r_dom*sin(theta_dom)*cos(phi_dom);
    y_cm_dom=r_dom*sin(theta_dom)*sin(phi_dom);
    z_cm_dom=r_dom*cos(theta_dom);

    x_temp2=rot_z[0]*x_cm_dom+rot_z[1]*y_cm_dom+rot_z[2]*z_cm_dom;
    y_temp2=rot_z[3]*x_cm_dom+rot_z[4]*y_cm_dom+rot_z[5]*z_cm_dom;
    z_temp2=rot_z[6]*x_cm_dom+rot_z[7]*y_cm_dom+rot_z[8]*z_cm_dom;

    x_temp=rot_y[0]*x_temp2+rot_y[1]*y_temp2+rot_y[2]*z_temp2;
    y_temp=rot_y[3]*x_temp2+rot_y[4]*y_temp2+rot_y[5]*z_temp2;
    z_temp=rot_y[6]*x_temp2+rot_y[7]*y_temp2+rot_y[8]*z_temp2;

    x_cm=rot_x[0]*x_temp+rot_x[1]*y_temp+rot_x[2]*z_temp;
    y_cm=rot_x[3]*x_temp+rot_x[4]*y_temp+rot_x[5]*z_temp;
    z_cm=rot_x[6]*x_temp+rot_x[7]*y_temp+rot_x[8]*z_temp;

    x_cm = x_cm - shift;
    
    theta = 0.0;
    error = 1000.0;

    if (y_cm >= 0.0) {
        theta_1 = 0.0;
        theta_2 = theta_HW;
    }
    if (y_cm < 0.0) {
        theta_1 = -theta_HW;
        theta_2 = 0.0;
    }

    while (error >= 0.01*d_theta){

        theta_rad = theta_1*c_PI/180.;
        theta_HW_rad = theta_HW*c_PI/180.;
        arg = c_PI*theta_rad/(2.0*theta_HW_rad);
        rx = R_t*(pow(abs(cos(arg)),n))*cos(theta_rad);
        ry = R_t*(pow(abs(cos(arg)),n))*sin(theta_rad);
        temp1 = R_t*cos(theta_rad)*n*(pow(abs(cos(arg)),n-1.))*sin(arg)*c_PI/(2.0*theta_HW_rad);
        temp2 = R_t*sin(theta_rad)*n*(pow(abs(cos(arg)),n-1.))*sin(arg)*c_PI/(2.0*theta_HW_rad);
        f_theta_1 = 2.0*((rx)-(x_cm))*(-(ry)-(temp1))+2.0*((ry)-(y_cm))*((rx)-(temp2));
        

        theta_rad = theta_2*c_PI/180.;
        theta_HW_rad = theta_HW*c_PI/180.;
        arg = c_PI*theta_rad/(2.0*theta_HW_rad);
        rx = R_t*(pow(abs(cos(arg)),n))*cos(theta_rad);
        ry = R_t*(pow(abs(cos(arg)),n))*sin(theta_rad);
        temp1 = R_t*cos(theta_rad)*n*(pow(abs(cos(arg)),n-1.))*sin(arg)*c_PI/(2.0*theta_HW_rad);
        temp2 = R_t*sin(theta_rad)*n*(pow(abs(cos(arg)),n-1.))*sin(arg)*c_PI/(2.0*theta_HW_rad);
        f_theta_2 = 2.0*((rx)-(x_cm))*(-(ry)-(temp1))+2.0*((ry)-(y_cm))*((rx)-(temp2));
        
        
        theta_center = (theta_1+theta_2)/2.0;
        
        theta_rad = theta_center*c_PI/180.;
        theta_HW_rad = theta_HW*c_PI/180.;
        arg = c_PI*theta_rad/(2.0*theta_HW_rad);
        rx = R_t*(pow(abs(cos(arg)),n))*cos(theta_rad);
        ry = R_t*(pow(abs(cos(arg)),n))*sin(theta_rad);
        temp1 = R_t*cos(theta_rad)*n*(pow(abs(cos(arg)),n-1.))*sin(arg)*c_PI/(2.0*theta_HW_rad);
        temp2 = R_t*sin(theta_rad)*n*(pow(abs(cos(arg)),n-1.))*sin(arg)*c_PI/(2.0*theta_HW_rad);
        f_theta_center = 2.0*((rx)-(x_cm))*(-(ry)-(temp1))+2.0*((ry)-(y_cm))*((rx)-(temp2));
        
        
        error = abs(theta_1-theta_center);
        if (f_theta_1*f_theta_center >= 0.0) {
            theta_1 = theta_center;
        }		
        if (f_theta_2*f_theta_center >= 0.0) {
            theta_2 = theta_center;
        }
    }


    theta = theta_center;
    if (abs(y_cm-0.0) <= epsilon1){
        theta = 0.0;
    }

    theta_rad = theta*c_PI/180.0;
    r_1 = R_t*(pow(abs(cos(c_PI*theta_rad/(2.0*theta_HW_rad))),n));
    r_2 = R_t*(pow(abs(cos(c_PI*(theta_rad+d_theta_rad)/(2.0*theta_HW_rad))),n));
    R_big = R_p*(pow(abs(cos(c_PI*theta_rad/(2.0*theta_HW_rad))),n));
    rx = r_1*cos(theta_rad);
    ry = r_1*sin(theta_rad);
    D = sqrt((rx-x_cm)*(rx-x_cm)+(ry-y_cm)*(ry-y_cm)+(z_cm)*(z_cm));
    if (D<=R_big){	   
        n_vec_i = r_1*cos(theta_rad)-r_2*cos(theta_rad+d_theta_rad);
        n_vec_j = r_1*sin(theta_rad)-r_2*sin(theta_rad+d_theta_rad);
        n_vec_mag = sqrt(n_vec_i*n_vec_i+n_vec_j*n_vec_j);
        n_hat_i = n_vec_i/n_vec_mag;
        n_hat_j = n_vec_j/n_vec_mag;
        r_dash = D;
        x_dash = 0.0;
        z_dash = z_cm;
        
        if (abs(y_cm/x_cm) <= abs(ry/rx)){
            y_dash = -sqrt(abs(r_dash*r_dash-z_dash*z_dash));
        } else {
            y_dash = sqrt(abs(r_dash*r_dash-z_dash*z_dash));
        }
        if (abs((y_cm-0.0)/c_SR) <= epsilon1) {
            if ((rx)*(rx) >= x_cm*x_cm){
                y_dash = -sqrt(abs(r_dash*r_dash-z_dash*z_dash));
            } else { 
                y_dash = sqrt(abs(r_dash*r_dash-z_dash*z_dash));	
            }
        }
        phi = atan2(z_dash,y_dash);
        B0_local = B0*log(T0*T0*R_p*R_p + 1.0)/log(T0*T0*R_big*R_big + 1.0);
        // If radius of curvature needed
        // double a_axis = M_PI/2.0/theta_HW_rad;
        // double r_axis = r_1;
        // double r_dash_axis = -n*a_axis*R_t*pow(cos(a_axis*theta_rad),n-1.)*sin(a_axis*theta_rad);
        // double r_double_dash_axis = a_axis*a_axis*n*R_t*(pow(cos(a_axis*theta_rad),n) - (n-1.)*pow(cos(a_axis*theta_rad),n-2.)*pow(sin(a_axis*theta_rad),2.0));
        // R0_Vandas = pow(r_axis*r_axis + r_dash_axis*r_dash_axis,1.5)/(r_axis*r_axis + 2.0*r_dash_axis*r_dash_axis - r_axis*r_double_dash_axis);

        // Another ad-hoc approximation for the radius of curvature
        R0_Vandas = (1.0+(theta_HW-45.0)/100)*R_p*(1.7-1.4*(R_p/R_t));
        
        //Gold-Hoyle (Only valid in cylindrical geometry)
        //B_x_dash = B0_local/(1.0+T0*T0*r_dash*r_dash);
        //B_y_dash = -ParseInputs::get_CME_helicity_sign()*B0_local*T0*r_dash*sin(phi)/(1.0+T0*T0*r_dash*r_dash);
        //B_z_dash =  ParseInputs::get_CME_helicity_sign()*B0_local*T0*r_dash*cos(phi)/(1.0+T0*T0*r_dash*r_dash);
        
        //Uniform twist solution of Vandas et al 2017
        B_x_dash = B0_local/(1.0+T0*T0*r_dash*r_dash);
        B_y_dash = -ParseInputs::get_CME_helicity_sign()*B0_local*R0_Vandas*T0*r_dash*sin(phi)/(1.0+T0*T0*r_dash*r_dash)/(R0_Vandas+r_dash*cos(phi));
        B_z_dash =  ParseInputs::get_CME_helicity_sign()*B0_local*R0_Vandas*T0*r_dash*cos(phi)/(1.0+T0*T0*r_dash*r_dash)/(R0_Vandas+r_dash*cos(phi));
        if ((R0_Vandas+r_dash*cos(phi)) <= 0.0) {
            B_x_dash = 0.0;
            B_y_dash = 0.0;
            B_z_dash = 0.0;
        }
        //Zeroth - order solution of Vandas et al 2017
        //B_x_dash = B0_local*R0_Vandas/(1.0+T0*T0*r_dash*r_dash)/(R0_Vandas+r_dash*cos(phi));
        //B_y_dash = -ParseInputs::get_CME_helicity_sign()*B0_local*R0_Vandas*T0*r_dash*sin(phi)/(1.0+T0*T0*r_dash*r_dash)/(R0_Vandas+r_dash*cos(phi));
        //B_z_dash =  ParseInputs::get_CME_helicity_sign()*B0_local*R0_Vandas*T0*r_dash*cos(phi)/(1.0+T0*T0*r_dash*r_dash)/(R0_Vandas+r_dash*cos(phi));
        
        
        
        if ((n_hat_j) >=0){
            t = acos(n_hat_i);
        } else { 
            t = -acos(n_hat_i);
        }

        B_x = B_x_dash*cos(t)-B_y_dash*sin(t);
        B_y = B_y_dash*cos(t)+B_x_dash*sin(t);
        B_z = B_z_dash;
        
        x_temp=rot_x_inv[0]*B_x+rot_x_inv[1]*B_y+rot_x_inv[2]*B_z;
        y_temp=rot_x_inv[3]*B_x+rot_x_inv[4]*B_y+rot_x_inv[5]*B_z;
        z_temp=rot_x_inv[6]*B_x+rot_x_inv[7]*B_y+rot_x_inv[8]*B_z;
        
        x_temp2=rot_y_inv[0]*x_temp+rot_y_inv[1]*y_temp+rot_y_inv[2]*z_temp;
        y_temp2=rot_y_inv[3]*x_temp+rot_y_inv[4]*y_temp+rot_y_inv[5]*z_temp;
        z_temp2=rot_y_inv[6]*x_temp+rot_y_inv[7]*y_temp+rot_y_inv[8]*z_temp;           

        B_x_dom=rot_z_inv[0]*x_temp2+rot_z_inv[1]*y_temp2+rot_z_inv[2]*z_temp2;
        B_y_dom=rot_z_inv[3]*x_temp2+rot_z_inv[4]*y_temp2+rot_z_inv[5]*z_temp2;
        B_z_dom=rot_z_inv[6]*x_temp2+rot_z_inv[7]*y_temp2+rot_z_inv[8]*z_temp2;
        
        sqBTotal = d_1_4PI*(B_x_dom*B_x_dom+B_y_dom*B_y_dom+B_z_dom*B_z_dom);
        e0TD = half*sqBTotal;
        
        //If an expanding velocity is needed
        a_hat_i = (x_cm-rx)/D;
        a_hat_j = (y_cm-ry)/D;
        a_hat_k = z_cm/D;
        D1 = sqrt((x_cm)*(x_cm)+(y_cm)*(y_cm)+(z_cm)*(z_cm));
        b_hat_i = x_cm/D1;
        b_hat_j = y_cm/D1;
        b_hat_k = z_cm/D1;
        V_rad = ParseInputs::get_CME_apex_speed()/(1.0+(R_p/R_t));
        
        V_exp = (D/R_t)*V_rad;
        V_x = V_rad*b_hat_i + V_exp*a_hat_i;
        V_y = V_rad*b_hat_j + V_exp*a_hat_j;
        V_z = V_rad*b_hat_k + V_exp*a_hat_k;
        
        x_temp=rot_x_inv[0]*V_x+rot_x_inv[1]*V_y+rot_x_inv[2]*V_z;
        y_temp=rot_x_inv[3]*V_x+rot_x_inv[4]*V_y+rot_x_inv[5]*V_z;
        z_temp=rot_x_inv[6]*V_x+rot_x_inv[7]*V_y+rot_x_inv[8]*V_z;
        
        x_temp2=rot_y_inv[0]*x_temp+rot_y_inv[1]*y_temp+rot_y_inv[2]*z_temp;
        y_temp2=rot_y_inv[3]*x_temp+rot_y_inv[4]*y_temp+rot_y_inv[5]*z_temp;
        z_temp2=rot_y_inv[6]*x_temp+rot_y_inv[7]*y_temp+rot_y_inv[8]*z_temp;            

        V_x_dom=rot_z_inv[0]*x_temp2+rot_z_inv[1]*y_temp2+rot_z_inv[2]*z_temp2;
        V_y_dom=rot_z_inv[3]*x_temp2+rot_z_inv[4]*y_temp2+rot_z_inv[5]*z_temp2;
        V_z_dom=rot_z_inv[6]*x_temp2+rot_z_inv[7]*y_temp2+rot_z_inv[8]*z_temp2;
        
        //If a simple radial velocity is needed
        // V_x_dom=ParseInputs::get_CME_apex_speed()*x_cm_dom/sqrt(x_cm_dom*x_cm_dom + y_cm_dom*y_cm_dom + z_cm_dom*z_cm_dom);
        // V_y_dom=ParseInputs::get_CME_apex_speed()*y_cm_dom/sqrt(x_cm_dom*x_cm_dom + y_cm_dom*y_cm_dom + z_cm_dom*z_cm_dom);
        // V_z_dom=ParseInputs::get_CME_apex_speed()*z_cm_dom/sqrt(x_cm_dom*x_cm_dom + y_cm_dom*y_cm_dom + z_cm_dom*z_cm_dom);

        a_CME(0) = ParseInputs::get_CME_density();
        a_CME(1) = V_x_dom*1.0e5;
        a_CME(2) = V_y_dom*1.0e5;
        a_CME(3) = V_z_dom*1.0e5;
        a_CME(4) = 0.0;//ParseInputs::get_CME_energy_control()*e0TD;
        a_CME(5) = B_x_dom;
        a_CME(6) = B_y_dom;
        a_CME(7) = B_z_dom;        
    }
        
}
PROTO_KERNEL_END(define_CMEF, define_CME)

template<
    typename T,
    template<MemType> typename MAP,
    MemType MEM>
inline void
BoxOp_EulerCubedSphere<T,MAP,MEM>::
Insert_CME(MBLevelBoxData<T, NUMCOMPS, MEM> &a_JU,
    MBLevelBoxData<T, 1, MEM> &a_dVolrLev,
    MBInterpOp &a_iop,
    MBDisjointBoxLayout &a_layout,
    double a_time,
    double a_dt,
    double a_gamma)
{
    if (ParseInputs::get_CME_type() == 1){
        if (procID() == 0) cout << "Inserting CME" << endl;
        HDF5Handler h5;
        auto map = CubedSphereShell::Map(a_JU.layout(), ghost());
        auto eulerOp = CubedSphereShell::Operator<BoxOp_EulerCubedSphere, double, HOST>(map);
        MBLevelBoxData<double, NUMCOMPS, HOST> USph_CME(a_layout, ghost());
        MBLevelBoxData<double, NUMCOMPS, HOST> Wout(a_layout, Point::Zeros());
        a_JU.exchange();
        a_JU.copyTo(USph_CME);
        CubedSphereShell::consToSphInterpEuler(USph_CME,a_iop,a_dVolrLev, 4); // Transform a_JU -> USph.
        int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
        int thetaCoord = (rCoord + 1) % DIM;
        int phiCoord = (rCoord + 2) % DIM;
        MBLevelBoxData<double, 8, HOST> dstData(a_layout, Point::Basis(rCoord) + NGHOST*Point::Basis(thetaCoord) + NGHOST*Point::Basis(phiCoord));
        if (ParseInputs::get_init_condition_type() == 3) BC_global.BoxData_to_BC<double,HOST>(dstData, map, a_time);
        for (auto dit : a_layout)
        {
            auto block = a_layout.block(dit);
            Box blockBox = a_layout.getBlock(block).domain().box();
            eulerOp[dit].PreStagePatch(USph_CME[dit],dstData[dit],a_dVolrLev[dit],blockBox,a_time,a_dt,0);
            Array<double, DIM> dx = eulerOp[dit].dx();
            double dxiPerp = dx[2];
            BoxData<double> radius(a_dVolrLev[dit].box());
            BoxData<double, DIM, HOST> Dr(a_dVolrLev[dit].box());
            BoxData<double, DIM, HOST> adjDr(a_dVolrLev[dit].box());
            eulerOp[dit].radialMetrics(radius, Dr, adjDr, a_dVolrLev[dit], Dr.box());
            
            BoxData<double, NUMCOMPS, HOST> W(USph_CME[dit].box());
            BoxData<double, NUMCOMPS, HOST> WBar(USph_CME[dit].box());
            eulerOp[dit].sphToPrim(W,WBar,USph_CME[dit],block); // Transform to primitive variables, spherical velocities.
            BoxData<double, NUMCOMPS, HOST> W_point(USph_CME[dit].box());
            W.copyTo(W_point);
            W_point += Stencil<double>::Laplacian()(W,-1./24.);  
            BoxData<double, NUMCOMPS, HOST> W_CME(USph_CME[dit].box());
            W_CME.setVal(0.);
            double half = 0.5;
            BoxData<double, DIM, HOST> XCart = forall_p<double,DIM,HOST>
                (f_cubedSphereMap3,radius.box(),radius,dx,half,half,block);  
            forallInPlace_p(define_CME,W_CME,XCart); // Set CME values in cartesian coordinates
            // Calculate A_matrix
            Box bx = W_CME.box();
            
            Array<Array<uint,DIM>,6> permute = {{1,2,0},{1,2,0},{1,0,2},{0,1,2},{1,0,2},{0,1,2}};
            Array<Array<int,DIM>,6> sign = {{1,-1,-1},{1,1,1},{1,-1,1},{1,1,1},{-1,1,1},{-1,-1,1}}; 
            Point high = bx.high();
            Point low = bx.low();
            high[0] = low[0];
            Box bx0(low,high);
            BoxData<double ,DIM,HOST,DIM> A_matrix(bx);
            double offseta = half;
            double offsetb = half;
            forallInPlace_p(f_Amatrix,bx0,A_matrix,permute[block],sign[block],
                            dxiPerp,offseta,offsetb);
            spreadSlice(A_matrix);

            BoxData<double,DIM,HOST,DIM> A_matrix_inv(bx);
            A_matrix_inv.setToZero();
            forallInPlace(f_matinv3by3,A_matrix_inv,A_matrix);

            BoxData<double,DIM,HOST> V_CME_sph(bx);
            BoxData<double,DIM,HOST> B_CME_sph(bx);

            BoxData<double,DIM,HOST> V_CME_cart = slice<double,NUMCOMPS,DIM,HOST>(W_CME,iVX);
            BoxData<double,DIM,HOST> B_CME_cart = slice<double,NUMCOMPS,DIM,HOST>(W_CME,iBX);
            double one = 1.0;  
            forallInPlace(f_matVecProd,V_CME_sph,A_matrix_inv,V_CME_cart,one);
            forallInPlace(f_matVecProd,B_CME_sph,A_matrix_inv,B_CME_cart,one);

            forallInPlace([ ] PROTO_LAMBDA
                        (Var<double, NUMCOMPS, HOST> &a_a_W,
                        Var<double, DIM, HOST> &a_V_sph,
                        Var<double, DIM, HOST> &a_B_sph)
            {  
                a_a_W(iVX) = a_V_sph(0);
                a_a_W(iVY) = a_V_sph(1);
                a_a_W(iVZ) = a_V_sph(2);
                a_a_W(iBX) = a_B_sph(0);
                a_a_W(iBY) = a_B_sph(1);
                a_a_W(iBZ) = a_B_sph(2);
            },W_CME, V_CME_sph, B_CME_sph);
            W_point.copyTo(Wout[dit]);
            // if (procID()==5) h5.writePatch(1,Wout[dit],"W_CME1");
            W_CME.copyTo(Wout[dit]);
            // if (procID()==5) h5.writePatch(1,Wout[dit],"W_CME2");
            forallInPlace([ ] PROTO_LAMBDA
                        (Var<double, NUMCOMPS, HOST> &a_W,
                        Var<double, NUMCOMPS, HOST> &a_W_CME,
                        double a_gamma)
            { 
                if (a_W_CME(iRHO) != 0){
                // a_W(iRHO) += a_W_CME(iRHO); 
                a_W(iVX) = a_W_CME(iVX);
                a_W(iVY) = a_W_CME(iVY);
                a_W(iVZ) = a_W_CME(iVZ);
                a_W(iBX) = a_W_CME(iBX);
                a_W(iBY) = a_W_CME(iBY);
                a_W(iBZ) = a_W_CME(iBZ);
                // a_W(iP) = 0;
                }
            },W_point, W_CME,a_gamma);
            MPI_Barrier(MPI_COMM_WORLD);

            W_point.copyTo(Wout[dit]);
            if (procID()==5) h5.writePatch(1,Wout[dit],"W_CME3");

            BoxData<double, NUMCOMPS, HOST> JUTemp;
            eulerOp[dit].primToCons(JUTemp, W_point, a_dVolrLev[dit], a_gamma, dx[2], block);
            if (procID()==5) h5.writePatch(1,JUTemp,"JUTemp");
            JUTemp.copyTo(a_JU[dit]);

        }
        a_JU.copyTo(USph_CME);
        CubedSphereShell::consToSphInterpEuler(USph_CME,a_iop,a_dVolrLev, 4); // Transform a_JU -> USph.
        for (auto dit : a_layout)
        {
            auto block = a_layout.block(dit);
            Box blockBox = a_layout.getBlock(block).domain().box();
            eulerOp[dit].PreStagePatch(USph_CME[dit],dstData[dit],a_dVolrLev[dit],blockBox,a_time,a_dt,0);
            Array<double, DIM> dx = eulerOp[dit].dx();
            double dxiPerp = dx[2];
            BoxData<double> radius(a_dVolrLev[dit].box());
            BoxData<double, DIM, HOST> Dr(a_dVolrLev[dit].box());
            BoxData<double, DIM, HOST> adjDr(a_dVolrLev[dit].box());
            eulerOp[dit].radialMetrics(radius, Dr, adjDr, a_dVolrLev[dit], Dr.box());
            
            BoxData<double, NUMCOMPS, HOST> W_point(USph_CME[dit].box());
            eulerOp[dit].sphToPrim_2O(W_point,USph_CME[dit],block, true); // Transform to primitive variables, spherical velocities.
            W_point.copyTo(Wout[dit]);
            if (procID()==5) h5.writePatch(1,Wout[dit],"W_CME4");
        }
      }
}