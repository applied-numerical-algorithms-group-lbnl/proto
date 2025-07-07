import h5py
import numpy as np
from datetime import datetime as dt
import time
from glob import glob
import os

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

#Initialize
r0 = 0.1
num_components = 8
time = toYearFraction(dt(2024, 3, 16, 12, 0, 0))
phys_domain = [0.1, 0, 0, 0.1, 6.28319, 3.14159]
Out_file = "swig_wsa2_smooth5_r21_5_bc_single.h5"
step_const = [0, 1, 0]
phi_face = True

# Get the directory containing the current Python file
current_directory = os.path.dirname(os.path.abspath(__file__))

# Move one level up
parent_directory = os.path.abspath(os.path.join(current_directory, '..'))

# Construct the new relative path
address_of_files = os.path.join(parent_directory, 'swig_wsa2_smooth5_r21_5_bc_single/')

files = glob(address_of_files+'rho_r1_idx*.h5')
files=sorted(files)
rho_filename = address_of_files+"rho_r1_idx.h5"
V_filename = address_of_files+"vr_r1_idx.h5"
T_filename = address_of_files+"t_r1_idx.h5"
B_filename = address_of_files+"br_r1_idx.h5"

#If n time dependent BCs availaible, they should be numbered like rho_upper_corona_Np1, ..., rho_upper_corona_Npn
print("Number of files:",len(files))
for i in range(len(files)):
    print("Converting file:",i+1)
    if (len(files) > 1):
        #make i+1 a string with 6 digits. Add leading zeros if necessary
        i_str = str(i+1813)
        i_str = i_str.zfill(6)
        rho_filename = address_of_files+"rho_r1_idx"+i_str+".h5"
        V_filename = address_of_files+"vr_r1_idx"+i_str+".h5"
        T_filename = address_of_files+"t_r1_idx"+i_str+".h5"
        B_filename = address_of_files+"br_r1_idx"+i_str+".h5"
    with h5py.File(rho_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        data_rho = list(f[a_group_key])
        data_rho = np.array(data_rho)
        # Step 1: Average adjacent elements to reduce dimensions from (181, 361) to (180, 360)
        # First, handle rows
        data_rho_reduced_rows = (data_rho[:-1, :] + data_rho[1:, :]) / 2.0
        # Since data_rho_reduced_rows is now (180, 361), we need to average adjacent columns
        if phi_face:
            # dont average, simply ignore the last column
            data_rho_reduced = data_rho_reduced_rows[:, :-1]
        else:
            data_rho_reduced = (data_rho_reduced_rows[:, :-1] + data_rho_reduced_rows[:, 1:]) / 2.0
        # Step 2: Transpose the array to get the final shape (360, 180)
        data_rho = data_rho_reduced.T
        dim_siz1 = data_rho.shape[0]
        dim_siz2 = data_rho.shape[1]
        data_rho = data_rho.flatten('F')
        b_group_key = list(f.keys())[2]
        theta = list(f[b_group_key])
        theta = np.array(theta)
        theta_reduced = (theta[:-1] + theta[1:]) / 2.0
        theta = theta_reduced.flatten('F')
        c_group_key = list(f.keys())[1]
        phi = list(f[c_group_key])
        phi = np.array(phi)
        if phi_face:
            # dont average, simply ignore the last column
            phi_reduced = phi[:-1]
        else:
            phi_reduced = (phi[:-1] + phi[1:]) / 2.0
        phi = phi_reduced.flatten('F')

    with h5py.File(V_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        print(list(f.keys()))
        data_Vr = list(f[a_group_key])
        data_Vr = np.array(data_Vr)
        # Step 1: Average adjacent elements to reduce dimensions from (181, 361) to (180, 360)
        # First, handle rows
        data_Vr_reduced_rows = (data_Vr[:-1, :] + data_Vr[1:, :]) / 2.0
        # Since data_Vr_reduced_rows is now (180, 361), we need to average adjacent columns
        if phi_face:
            # dont average, simply ignore the last column
            data_Vr_reduced = data_Vr_reduced_rows[:, :-1]
        else:
            data_Vr_reduced = (data_Vr_reduced_rows[:, :-1] + data_Vr_reduced_rows[:, 1:]) / 2.0
        # Step 2: Transpose the array to get the final shape (360, 180)
        data_Vr = data_Vr_reduced.T
        data_Vr = data_Vr.flatten('F')
        data_Vr = data_Vr*1e5  #make it cm/s
        data_Vp = data_Vr*0.0
        data_Vt = data_Vr*0.0

    with h5py.File(T_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        data_T = list(f[a_group_key])
        data_T = np.array(data_T)
        # Step 1: Average adjacent elements to reduce dimensions from (181, 361) to (180, 360)
        # First, handle rows
        data_T_reduced_rows = (data_T[:-1, :] + data_T[1:, :]) / 2.0
        # Since data_T_reduced_rows is now (180, 361), we need to average adjacent columns
        if phi_face:
            # dont average, simply ignore the last column
            data_T_reduced = data_T_reduced_rows[:, :-1]
        else:
            # Average adjacent columns
            data_T_reduced = (data_T_reduced_rows[:, :-1] + data_T_reduced_rows[:, 1:]) / 2.0
        # Step 2: Transpose the array to get the final shape (360, 180)
        data_T = data_T_reduced.T
        data_T = data_T.flatten('F')   

    data_P = 2.0*data_rho*1.3806505e-16*data_T/1e-12   # in pico Dyne

    with h5py.File(B_filename, "r") as f:
        a_group_key = list(f.keys())[0]
        data_Br = list(f[a_group_key])
        data_Br = np.array(data_Br)
        # Step 1: Average adjacent elements to reduce dimensions from (181, 361) to (180, 360)
        # First, handle rows
        data_Br_reduced_rows = (data_Br[:-1, :] + data_Br[1:, :]) / 2.0
        # Since data_Br_reduced_rows is now (180, 361), we need to average adjacent columns
        if phi_face:
            # dont average, simply ignore the last column
            data_Br_reduced = data_Br_reduced_rows[:, :-1]
        else:
            # Average adjacent columns
            data_Br_reduced = (data_Br_reduced_rows[:, :-1] + data_Br_reduced_rows[:, 1:]) / 2.0
        # Step 2: Transpose the array to get the final shape (360, 180)
        data_Br = data_Br_reduced.T
        data_Br = data_Br.flatten('F')
        data_Br = data_Br*1e6  #make it micro Gauss
        data_Bp = data_Br*0.0
        data_Bt = data_Br*0.0 

    data_final = np.concatenate((data_rho,data_Vr,data_Vp, data_Vt, data_P, data_Br, data_Bp, data_Bt))
    dtheta = np.full(dim_siz2, np.pi/dim_siz2)
    if (i==0):
        with h5py.File(Out_file, "w") as data_file:
            data_file.create_dataset("data"+str(i), data=data_final,dtype='float64')
            data_file["data"+str(i)].attrs["time"] = i
            data_file.attrs['domain'] = [dim_siz1,dim_siz2]
            data_file.attrs['num_components'] = num_components
            data_file.attrs['num_datasets'] = len(files)
            data_file.attrs['time'] = time
            data_file.attrs['r0'] = r0
            data_file.attrs['datasets_time'] = range(len(files))
            grp = data_file.create_group("geometry")
            grp.attrs['phys_domain'] = phys_domain
            grp.attrs['step_const'] = step_const
            grp.create_dataset("dtheta", data=dtheta,dtype='float64')
            grp.create_dataset("theta", data=theta,dtype='float64')
            grp.create_dataset("phi", data=phi,dtype='float64')
    if (i!=0):
        with h5py.File(Out_file, "a") as data_file:
            data_file.create_dataset("data"+str(i), data=data_final,dtype='float64')
            data_file["data"+str(i)].attrs["time"] = i

print("Time in fraction of year:", time)
print("Conversion completed. Output file:", Out_file)
print("Output file saved in:", os.path.abspath(Out_file))
print("Data shape:", data_final.shape)
print("Data domain:", [dim_siz1, dim_siz2])
print("Data type:", data_final.dtype)
print("Number of components:", num_components)
print("Physical domain:", phys_domain)