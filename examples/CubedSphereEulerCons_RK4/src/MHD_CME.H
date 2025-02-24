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
        a_CME(4) = ParseInputs::get_CME_energy_control()*e0TD;
        a_CME(5) = B_x_dom;
        a_CME(6) = B_y_dom;
        a_CME(7) = B_z_dom;        
    }
        
}
PROTO_KERNEL_END(define_CMEF, define_CME)



template <typename T, MemType MEM>
PROTO_KERNEL_START
void superimpose_CMEF(const Point& a_pt,
                      Var<T, NUMCOMPS, MEM> &a_U,
                      Var<T, NUMCOMPS, MEM> &a_CME)
{
    if (abs(a_CME(1)) != 0){
        a_U(0) = a_U(0) + a_CME(0);			
        a_U(1) = a_U(0)*a_CME(1);			
        a_U(2) = a_U(0)*a_CME(2);			
        a_U(3) = a_U(0)*a_CME(3);	
        a_U(4) = a_U(4) + a_CME(4);	
        a_U(5) = a_CME(5);			
        a_U(6) = a_CME(6);			
        a_U(7) = a_CME(7);		
    }
}
PROTO_KERNEL_END(superimpose_CMEF, superimpose_CME)


void define_CME_calc(State& a_CME,
                     const Vector& a_x)
{
    forallInPlace_p(define_CME,a_CME,a_x);
}

// void Insert_CME(MHDLevelDataState& a_state,
//             const int a_k,
//             const double a_time,
//             const double a_dt)
// {
//     double physical_time = getPhysTime(a_time);

//     double a_dx = a_state.m_dx;
//     double a_dy = a_state.m_dy;
//     double a_dz = a_state.m_dz;
//     double a_gamma = a_state.m_gamma;

//     if (!a_state.m_CME_inserted && (physical_time >= ParseInputs::get_CME_Enter_Time())) {
//         if (procID() == 0) cout << "Inserting CME" << endl;
//         for (auto dit : a_state.m_CME){	
//             a_state.m_CME[dit].setVal(0.0);	  
//             Box dbx1 = a_state.m_CME[dit].box();
//             BoxData<double,DIM> x(dbx1);		
//             MHD_Mapping::get_sph_coords_cc(x, dbx1, a_dx, a_dy, a_dz);
//             forallInPlace_p(define_CME,a_state.m_CME[dit],x);
//             forallInPlace_p(superimpose_CME,a_state.m_U[dit],a_state.m_CME[dit]);
//         }
//         a_state.m_CME_inserted = true;
//         MHD_Output_Writer::Write_data(a_state, a_k-1, physical_time, a_dt, true);
//     }
// }

