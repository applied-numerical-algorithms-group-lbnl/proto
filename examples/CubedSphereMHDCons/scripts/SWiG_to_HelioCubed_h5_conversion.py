import h5py
import numpy as np
from datetime import datetime as dt
import time
from glob import glob
import os
import matplotlib.pyplot as plt

def toYearFraction(date):
    def sinceEpoch(date): return time.mktime(date.timetuple())
    s = sinceEpoch
    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)
    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    return date.year + yearElapsed/yearDuration

# --- Parker spiral helper ---
def R_omega_cms(r0_au, rotation_period_days= 25.38):
    """
    Returns -R*Omega (cm/s) for a given inner-boundary radius r0 in AU.
    Uses Carrington sidereal rotation period by default (25.38 days).
    """
    AU_cm = 1.495978707e13        # 1 AU IN cm
    R_cm = r0_au * AU_cm
    omega = 2.0*np.pi / (rotation_period_days * 86400.0)   # rad/s
    return R_cm * omega         # cm/s 

# Settings
save_png = True  # Toggle this to enable/disable PNG output
png_output_dir = "frame_images_corotating"
if save_png:
    os.makedirs(png_output_dir, exist_ok=True)

#Initialize
r0 = 0.1
num_components = 8
time = toYearFraction(dt(2024, 3, 16, 12, 0, 0))
First_file = 1813
phys_domain = [0.1, 0, 0, 0.1, 6.28319, 3.14159]
# Out_file = "psi_eclipse24_swig_wsa2_smth5_bc_r21_5_corotating.h5"
Out_file = "psi_eclipse24_swig_wsa2_smth5_bc_r21_5_corotating_Bp_two_frames.h5"
step_const = [0, 1, 0]
phi_face = True
time_cadence = 1.0/24.

Br_scale_factor = 1.0 # this is 2.0, in Hegde et al,2025, otherwise keep it as 1.0 just keep the strength as it is.
redV = 0.0 #km/s # this is 75 km/s, in Hegde et al,2025, otherwise keep it as 0.0 just keep the speed as it is
minV = 0.0  #km/s # this is 200 km/s, in Hegde et al,2025, otherwise keep it as 0.0 just keep the minimum speed as it is
maxV = 700000.0 #km/s  # this is 700 km/s, in Hegde et al,2025, otherwise keep it as 0.0 just keep the minimum speed as it is

current_directory = os.path.dirname(os.path.abspath(__file__))
Out_file = os.path.join(current_directory, Out_file)
parent_directory = os.path.abspath(os.path.join(current_directory, '..'))
address_of_files = os.path.join(parent_directory, 'psi_eclipse24_swig_wsa2_smth5_bc_r21_5_corotating/')

files = sorted(glob(address_of_files+'rho_r1_idx*.h5'))
#keep only the first 1 files for testing
files = files[:2]
print("Number of files:", len(files))

for i in range(len(files)):
    print("Converting file:", i+1)
    # if len(files) > 1:
    i_str = str(i + First_file).zfill(6)
    rho_filename = address_of_files + f"rho_r1_idx{i_str}.h5"
    V_filename = address_of_files + f"vr_r1_idx{i_str}.h5"
    T_filename = address_of_files + f"t_r1_idx{i_str}.h5"
    B_filename = address_of_files + f"br_r1_idx{i_str}.h5"

    # Read and process rho
    with h5py.File(rho_filename, "r") as f:
        data_rho = np.array(f[list(f.keys())[0]])
        data_rho = ((data_rho[:-1, :] + data_rho[1:, :]) / 2.0)
        if phi_face:
            data_rho = data_rho[:, :-1]
        else:
            data_rho = (data_rho[:, :-1] + data_rho[:, 1:]) / 2.0
        data_rho = data_rho.T
        dim_siz1, dim_siz2 = data_rho.shape
        rho_flat = data_rho.flatten('F')
        
        # Handling angles
        theta = np.array(f[list(f.keys())[2]])
        theta = ((theta[:-1] + theta[1:]) / 2.0).flatten('F')
        phi = np.array(f[list(f.keys())[1]])
        phi = (phi[:-1] if phi_face else (phi[:-1] + phi[1:]) / 2.0).flatten('F')

    # Read and process velocity
    with h5py.File(V_filename, "r") as f:
        data_Vr = np.array(f[list(f.keys())[0]])
        data_Vr = ((data_Vr[:-1, :] + data_Vr[1:, :]) / 2.0)
        data_Vr = (data_Vr[:, :-1] if phi_face else (data_Vr[:, :-1] + data_Vr[:, 1:]) / 2.0)#.T
        data_Vr_kms = data_Vr.T
        
        # Reducing radial velocity and clipping
        data_Vr_kms = np.clip(data_Vr_kms - redV, minV, maxV)
        data_Vr = data_Vr_kms* 1e5
        
        #data_Vr *= 1e5
        Vr_flat = data_Vr.flatten('F')
        Vp_flat = Vr_flat * 0.0
        Vt_flat = Vr_flat * 0.0

    # Read and process temperature
    with h5py.File(T_filename, "r") as f:
        data_T = np.array(f[list(f.keys())[0]])
        data_T = ((data_T[:-1, :] + data_T[1:, :]) / 2.0)
        data_T = (data_T[:, :-1] if phi_face else (data_T[:, :-1] + data_T[:, 1:]) / 2.0).T
        T_flat = data_T.flatten('F')

    # Compute pressure
    P_flat = 2.0 * rho_flat * 1.3806505e-16 * T_flat / 1e-12

    # Read and process B field
    with h5py.File(B_filename, "r") as f:
        data_Br = np.array(f[list(f.keys())[0]])
        data_Br = ((data_Br[:-1, :] + data_Br[1:, :]) / 2.0)
        data_Br = (data_Br[:, :-1] if phi_face else (data_Br[:, :-1] + data_Br[:, 1:]) / 2.0).T
        #data_Br *= 1e6
        # apply scale factor
        data_Br *= 1e6*Br_scale_factor
        Br_flat = data_Br.flatten('F')
        
        # Bphi Calculation
        R_OMEGA_cms = R_omega_cms(r0_au=r0) 
        # Build 2D sin(theta) grid
        theta_2d = np.tile(theta, (dim_siz1, 1))            # (dim_siz1, dim_siz2)
        sin_theta_flat = np.sin(theta_2d).flatten('F')      # length = dim_siz1*dim_siz2
        
        Bp_flat = -(R_OMEGA_cms * sin_theta_flat * Br_flat) / Vr_flat

        # For plotting
        data_Bp = np.reshape(Bp_flat, (dim_siz1, dim_siz2), order='F')
    
        #Adjust Br to conserve |B|, then restore original Br sign (C++ behavior)
        Br_sign = np.sign(Br_flat)
        Br_adj_mag = np.sqrt(Br_flat**2 - Bp_flat**2)
        Br_flat = Br_sign * Br_adj_mag
        
        #Bp_flat = Br_flat * 0.0
        Bt_flat = Br_flat * 0.0

    data_final = np.concatenate((rho_flat, Vr_flat, Vp_flat, Vt_flat, P_flat, Br_flat, Bp_flat, Bt_flat))
    dtheta = np.full(dim_siz2, np.pi / dim_siz2)

    if i == 0:
        with h5py.File(Out_file, "w") as data_file:
            data_file.create_dataset("data"+str(i), data=data_final)
            data_file["data"+str(i)].attrs["time"] = i * time_cadence
            data_file.attrs['domain'] = [dim_siz1, dim_siz2]
            data_file.attrs['num_components'] = num_components
            data_file.attrs['num_datasets'] = len(files)
            data_file.attrs['time'] = time
            data_file.attrs['r0'] = r0

            # data_file.attrs['datasets_time'] = [j * time_cadence for j in range(len(files))]
            # Create datasets_time as a DATASET (not an attribute)
            times = np.arange(len(files), dtype=np.float64) * time_cadence
            data_file.create_dataset("datasets_time", data=times)

            grp = data_file.create_group("geometry")
            grp.attrs['phys_domain'] = phys_domain
            grp.attrs['step_const'] = step_const
            grp.create_dataset("dtheta", data=dtheta)
            grp.create_dataset("theta", data=theta)
            grp.create_dataset("phi", data=phi)
    else:
        with h5py.File(Out_file, "a") as data_file:
            data_file.create_dataset("data"+str(i), data=data_final)
            data_file["data"+str(i)].attrs["time"] = i * time_cadence

    # Save PNG image if enabled
    if save_png:
        fig, axs = plt.subplots(2, 2, figsize=(22, 10))

        pcm0 = axs[0, 0].pcolormesh(data_rho.T, shading='auto')
        axs[0, 0].set_title(f"Density Frame {i}")
        axs[0, 0].invert_yaxis()
        fig.colorbar(pcm0, ax=axs[0, 0])

        pcm1 = axs[0, 1].pcolormesh(data_Vr.T, shading='auto')
        axs[0, 1].set_title(f"Radial Velocity Frame {i}")
        axs[0, 1].invert_yaxis()
        fig.colorbar(pcm1, ax=axs[0, 1])

        pcm2 = axs[1, 0].pcolormesh(data_Br.T, shading='auto')
        axs[1, 0].set_title(f"Radial Magnetic Field Frame {i}")
        axs[1, 0].invert_yaxis()
        fig.colorbar(pcm2, ax=axs[1, 0])
        
        #To check Bphi just comment our Temperature panel and uncomment this
        pcm3 = axs[1, 1].pcolormesh(data_Bp.T, shading='auto')
        axs[1, 1].set_title(f"Azimuthal Magnetic Field Frame {i_str}")
        axs[1, 1].invert_yaxis()
        fig.colorbar(pcm3, ax=axs[1, 1])

        # pcm3 = axs[1, 1].pcolormesh(data_T.T, shading='auto')
        # axs[1, 1].set_title(f"Temperature Frame {i}")
        # axs[1, 1].invert_yaxis()
        # fig.colorbar(pcm3, ax=axs[1, 1])

        plt.tight_layout()
        plt.savefig(os.path.join(png_output_dir, f"frame_{i:04d}.png"))
        plt.close()

print("Time in fraction of year:", time)
print("Conversion completed. Output file:", Out_file)
print("Output file saved in:", os.path.abspath(Out_file))
print("Data shape:", data_final.shape)
print("Data domain:", [dim_siz1, dim_siz2])
print("Data type:", data_final.dtype)
print("Number of components:", num_components)
print("Physical domain:", phys_domain)
if save_png:
    print(f"PNG images saved to: {png_output_dir}")
