#pragma once
#include "Proto.H"
#include "hdf5.h"
#include "Inputs_Parsing.H"

using namespace Proto;

class MHDReader
{
public:
	inline MHDReader() {}

	std::vector<BoxData<double, 8, HOST>> BC_data;
	std::vector<BoxData<double, 8, HOST>> BC_data_extended;
	std::vector<double> dtheta;
	std::vector<double> theta;
	std::vector<double> phi;
	Box Box_BC;
        Array<double,NUMCOMPS> m_scalevars =
                  {RHO00,VEL00*RHO00,VEL00*RHO00,VEL00*RHO00,P00,B00,B00,B00};

	inline void file_to_BoxData_vec(std::string a_filename)
	{
		PR_TIMERS("MHDReader::file_to_BoxData_vec()");
		int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
		int thetaCoord = (rCoord + 1) % 3;
		int phiCoord   = (rCoord + 2) % 3;

		if (procID() == 0) cout << "Reading BC data from file: " << a_filename << endl;

		// Geometry can be read by all (small)
		readGeom(dtheta, a_filename, /*parallel_open=*/true);
		readTheta(theta, a_filename, /*parallel_open=*/true);
		readPhi(phi, a_filename, /*parallel_open=*/true);
		Box_BC = readDomain(a_filename, /*parallel_open=*/true);

		// Build rotated BT box (same on all ranks)
		Point LT = Box_BC.low();
		LT[thetaCoord] = Box_BC.low()[phiCoord];
		LT[phiCoord]   = Box_BC.low()[thetaCoord];
		Point HT = Box_BC.high();
		HT[thetaCoord] = Box_BC.high()[phiCoord];
		HT[phiCoord]   = Box_BC.high()[thetaCoord];
		Box BT(LT, HT);

		// Extended box (same on all ranks)
		Point high = BT.high();
		high[phiCoord] = 2 * (BT.high()[phiCoord] + 1) - 1;
		Point low  = BT.low();
		low[phiCoord]  = -1;
		Box bx2(low, high);

		// Only rank 0 reads & prepares data frames
		if (procID() == 0)
		{
			// readData(BC_data, a_filename);
			readData(BC_data, a_filename, /*parallel_open=*/false);  // root-only serial open
			int force_single_frame = ParseInputs::get_BC_force_single_frame();
			if (force_single_frame == 1) BC_data.resize(1);

			// rotate theta<->phi into BT
			auto t = Point::Basis(thetaCoord);
			auto p = Point::Basis(phiCoord);
			CoordPermutation R{{t, p}, {p, t}};

			BC_data_extended.clear();
			BC_data_extended.resize(BC_data.size());

			for (int i = 0; i < (int)BC_data.size(); ++i)
			{
				BC_data[i].rotate(BT, R);

				// define & init extended on root only
				BC_data_extended[i].define(bx2);
				// for (int comp = 0; comp < 8; ++comp) {
				// 	BC_data_extended[i].setVal(m_scalevars[comp], bx2, comp);
				// }
				// wrap copy: [BT] at offset + direct copy
				BC_data[i].copyTo(BC_data_extended[i], BT, Point::Basis(phiCoord) * (BT.high()[phiCoord] + 1));
				BC_data[i].copyTo(BC_data_extended[i]);
			}
			cout << "Geometric + BC frames prepared on rank 0." << endl;
		}
		barrier();
		// All ranks keep these shapes for later
		Box_BC = BT;
	}

	// helper: broadcast a BoxData buffer (contiguous doubles)
	inline void bcast_boxdata(BoxData<double,8,HOST>& buf, int root=0)
	{
	#ifdef PR_MPI
		// Broadcast in chunks if ever larger than INT_MAX elements (rare)
		const size_t nElems = (size_t)buf.size();  // (#cells * 8)
		double* ptr = buf.data();
		const size_t INT_MAX_ELEMS = (size_t)std::numeric_limits<int>::max();
		size_t off = 0;
		while (off < nElems) {
			size_t chunk = std::min(nElems - off, INT_MAX_ELEMS);
			MPI_Bcast(ptr + off, (int)chunk, MPI_DOUBLE, root, MPI_COMM_WORLD);
			off += chunk;
		}
	#else
		(void)buf; (void)root; // no-op in serial
	#endif
	}

	// Using this instead of Stencil-based fractional shift for BoxData increased speed significantly
	template <typename T, MemType MEM>
	inline void
	apply_frac_shift_phi(BoxData<T,NUMCOMPS,MEM>       & a_out,
						 BoxData<T,NUMCOMPS,MEM> & a_in,
						int                            a_phiCoord,
						double                         a_frac)
	{
		// IMPORTANT: alias shift sign is opposite of Stencil::Shift
		// To mimic Shift(-Basis(phiCoord)), we use +Basis(phiCoord) in alias:
		auto in_left = alias(a_in, Point::Basis(a_phiCoord) * 1);

		forallInPlace(
			[=] PROTO_LAMBDA(
				Var<T,NUMCOMPS,MEM>       & l_a_out,
				const Var<T,NUMCOMPS,MEM> & l_a_in,
				const Var<T,NUMCOMPS,MEM> & l_in_left)
			{
				for (int c = 0; c < NUMCOMPS; ++c)
				{
					l_a_out(c) =
						(T(1.0) - T(a_frac)) * l_a_in(c)
					+ T(a_frac)            * l_in_left(c);
				}
			},
			a_out, a_in, in_left);
	}

	template <typename T, MemType MEM>
	inline void BoxData_to_BC(MBLevelBoxData<T, 8, MEM> &a_dstData,
							const MBLevelMap<MBMap_CubedSphereShell, MEM> &a_map,
							T a_time)
	{
		PR_TIMERS("MHDReader::BoxData_to_BC");
		HDF5Handler h5;
		int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
		int thetaCoord = (rCoord + 1) % 3;
		int phiCoord = (rCoord + 2) % 3;

		// --- minimal change: define boxes without touching BC_data[...] on non-root
		Box bx1 = Box_BC;  // was: BC_data[0].box()
		Point eHigh = bx1.high(); eHigh[phiCoord] = 2*(bx1.high()[phiCoord] + 1) - 1;
		Point eLow  = bx1.low();  eLow [phiCoord] = -1;
		Box bx2(eLow, eHigh);     // was: BC_data_extended[0].box()
		BoxData<double, 8, HOST> data_0(bx1);
		BoxData<double, 8, HOST> BC_data_temp1(bx1);
		BoxData<double, 8, HOST> BC_data_temp2(bx1);

		Point high = bx1.high();
		int domainSizez = high[phiCoord] + 1;
		// --- minimal change: broadcast num_frames so all ranks take same branch
		int num_frames = 0;
		if (procID()==0) num_frames = (int)BC_data_extended.size();
#ifdef PR_MPI
		MPI_Bcast(&num_frames, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

		int low_frame_now = (int)( (double)a_time / (ParseInputs::get_BC_cadence() * 3600.0) );
		if (num_frames == 1 || low_frame_now >= num_frames - 1)
		{

			BoxData<double, 8, HOST> data_a(bx2);
			// Rotate the data
			double carr_rot_time = 25.38 * 24 * 60 * 60; // Seconds
			// We should use sidereal time for this. 25.38 days. That's the rotation time from a fixed location.
			// Carrington rotation time (27.2753 days) is from Earth's prespective.
			double angle_to_rotate = fmod(360.0 * (double)a_time / carr_rot_time, 360.0);

			// Can use this if no rotation is needed
			// angle_to_rotate = 0;

			int cells_to_rotate = angle_to_rotate / (360 / domainSizez);
			double needed_fraction = angle_to_rotate / (360 / domainSizez) - cells_to_rotate;
			cells_to_rotate = cells_to_rotate % domainSizez;
			cells_to_rotate = domainSizez - cells_to_rotate;
			if (procID()==0) {
				if (num_frames == 1) BC_data_extended[0].copyTo(data_a, bx2, Point::Basis(phiCoord) * (-cells_to_rotate));
				if (low_frame_now >= num_frames - 1) BC_data_extended[num_frames - 1].copyTo(data_a, bx2, Point::Basis(phiCoord) * (-cells_to_rotate));
			}
			bcast_boxdata(data_a);
			apply_frac_shift_phi<T, MEM>(BC_data_temp1, data_a, phiCoord, needed_fraction);
			BC_data_temp1.copyTo(data_0);
		}
		else
		{
			double h5_cadence = ParseInputs::get_BC_cadence();
			double time_temp = a_time;
			int low_frame = time_temp / (h5_cadence * 3600);
			int next_frame = low_frame + 1;
			PROTO_ASSERT(low_frame < num_frames,
									 "MHD_Set_Boundary_Values::interpolate_h5_BC | Error: next BC frame not available. %i >= %lu", low_frame, num_frames);
			double needed_fraction = (time_temp - (low_frame * h5_cadence * 3600)) / (h5_cadence * 3600);

			BoxData<double, 8, HOST> data_a(bx2);
			BoxData<double, 8, HOST> data_b(bx2);
                     
			// Rotate the data
			double carr_rot_time = 25.38 * 24 * 60 * 60; // Seconds
			// We should use sidereal time for this. 25.38 days. That's the rotation time from a fixed location.
			// Carrington rotation time (27.2753 days) is from Earth's prespective.
			
			int frames_rotate = ParseInputs::get_BC_frame_rotate();

			double effective_time1, effective_time2, angle_to_rotate1, angle_to_rotate2, needed_fraction1, needed_fraction2;
			int cells_to_rotate1, cells_to_rotate2;

			if (frames_rotate == 1){
				// Use this when frames rotate. Tae's WSA BCs have rotating frames.
				effective_time1 = 360 * (time_temp + low_frame * h5_cadence * 3600) / carr_rot_time;
				effective_time2 = 360 * (time_temp + next_frame * h5_cadence * 3600) / carr_rot_time; // Should effective time be the same for both? No, if frames rotate as well, Yes, if frames are fixed.
				angle_to_rotate1 = fmod(effective_time1, 360);
				angle_to_rotate2 = fmod(effective_time2, 360);
			}

			if (frames_rotate == 0){
				//In SWiG, frames are fixed. In that case, use this
				angle_to_rotate1 = fmod(360 * a_time / carr_rot_time, 360);
				angle_to_rotate2 = fmod(360 * a_time / carr_rot_time, 360);
			}

			cells_to_rotate1 = angle_to_rotate1 / (360 / domainSizez);
			cells_to_rotate2 = angle_to_rotate2 / (360 / domainSizez);
			needed_fraction1 = angle_to_rotate1 / (360 / domainSizez) - cells_to_rotate1;
			needed_fraction2 = angle_to_rotate2 / (360 / domainSizez) - cells_to_rotate2;
			cells_to_rotate1 = cells_to_rotate1 % domainSizez;
			cells_to_rotate2 = cells_to_rotate2 % domainSizez;
			cells_to_rotate1 = domainSizez - cells_to_rotate1;
			cells_to_rotate2 = domainSizez - cells_to_rotate2;

			// if(procID()==0) cout << " low_frame = " << low_frame << " needed_fraction = " << needed_fraction << " cells_to_rotate1 " << cells_to_rotate1 << " cells_to_rotate2 " << cells_to_rotate2 << endl;

			if (procID()==0) {
				PR_TIMERS("MHDReader::BoxData_to_BC | copyTo frames");
				BC_data_extended[low_frame].copyTo(data_a, bx2, Point::Basis(phiCoord) * (-cells_to_rotate1));
				BC_data_extended[next_frame].copyTo(data_b, bx2, Point::Basis(phiCoord) * (-cells_to_rotate2));
			}
			{
				PR_TIMERS("MHDReader::BoxData_to_BC | bcast frames");
			bcast_boxdata(data_a);
			bcast_boxdata(data_b);
			}

			{ PR_TIMERS("MHDReader::BoxData_to_BC | apply_frac_shift_phi");
			apply_frac_shift_phi<T, MEM>(BC_data_temp1, data_a, phiCoord, needed_fraction1);
			apply_frac_shift_phi<T, MEM>(BC_data_temp2, data_b, phiCoord, needed_fraction2);
			}
			{ PR_TIMERS("MHDReader::BoxData_to_BC | final interpolation");
			forallInPlace([] PROTO_LAMBDA(Var<T, NUMCOMPS, MEM> & a_BC_data_interpolated,
										Var<T, NUMCOMPS, MEM> & a_BC_data_lo,
										Var<T, NUMCOMPS, MEM> & a_BC_data_hi,
										double a_fraction)
			{                    
                for (int i=0; i<NUMCOMPS; i++){
                	a_BC_data_interpolated(i) = a_fraction*a_BC_data_hi(i) + (1.0-a_fraction)*a_BC_data_lo(i);
                } 
			},
				data_0, BC_data_temp1, BC_data_temp2, needed_fraction);
			}
		}

		// Making use of cached BC operator to avoid rebuild each time, much faster
		// --- Cached source layout / data / BC operator ---
		struct BCCache
		{
			bool initialized = false;
			MBDisjointBoxLayout srcLayout;
			std::unique_ptr<MBLevelBoxData<double, 8, HOST>> srcData;
			MBInterpOp BCOp;
			int cachedNumProc = -1;
		};

		static BCCache s_cache;

		if (!s_cache.initialized || s_cache.cachedNumProc != numProc())
		{
			// Rebuild if first time, or if numProc() changed (paranoia)
			MBProblemDomain srcDomain(numProc());
			for (BlockIndex bi = 0; bi < numProc(); bi++)
			{
				srcDomain.defineDomain(bi, Box_BC.sizes());
			}

			s_cache.srcLayout = MBDisjointBoxLayout(srcDomain, Box_BC.sizes());
			s_cache.srcData   = std::make_unique<MBLevelBoxData<double, 8, HOST>>(
									s_cache.srcLayout, Point::Zeros());

			// Build interpolation op once, based on the map + srcLayout geometry
			s_cache.BCOp = CubedSphereShell::BCBilinearOp(
				a_map, s_cache.srcLayout, dtheta, theta, phi, Side::Lo, NGHOST);

			s_cache.cachedNumProc = numProc();
			s_cache.initialized   = true;
		}

		// Copy BC "image" into the (one) source patch this operator expects
		data_0.copyTo((*s_cache.srcData)[*s_cache.srcLayout.begin()]);

		a_dstData.setVal(0);
		{
			PR_TIMERS("MHDReader::BoxData_to_BC::applyBCOp");
			s_cache.BCOp.apply(a_dstData, *s_cache.srcData);
		}
	}


	// add a bool flag to choose serial vs parallel open
	inline void readFile(hid_t* a_file, const std::string& a_filename, bool parallel_open)
	{
	#ifdef PR_MPI
		if (parallel_open) {
			// collective open on WORLD (call this ONLY when all ranks participate)
			hid_t p_access = H5Pcreate(H5P_FILE_ACCESS);
			assert(p_access >= 0);
			assert(H5Pset_fapl_mpio(p_access, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
			*a_file = H5Fopen(a_filename.c_str(), H5F_ACC_RDONLY, p_access);
			assert(H5Pclose(p_access) >= 0);
		} else {
			// serial open (safe to call on rank 0 only)
			*a_file = H5Fopen(a_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		}
	#else
		*a_file = H5Fopen(a_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	#endif
		if (*a_file < 0) {
			pr_out() << "H5Fopen failed for: " << a_filename << "\n";
	#ifdef PR_MPI
			MPI_Abort(MPI_COMM_WORLD, 1);
	#else
			std::abort();
	#endif
		}
	}

	// small helper to avoid world barriers inside root-only paths
	inline void barrier_world_if_all(bool do_barrier) {
	#ifdef PR_MPI
		if (do_barrier) MPI_Barrier(MPI_COMM_WORLD);
	#endif
	}

	inline Box _readDomain(hid_t *a_file)
	{
		// READ NUM COMPONENTS
		// int numComps;
		// auto att_comps = H5Aopen(*a_file, "num_components", H5P_DEFAULT);
		// assert(H5Aread(att_comps, H5T_STD_I32LE, &numComps) >= 0);
		// assert(H5Aclose(att_comps) >= 0);

		// PROTO_ASSERT(numComps == C,
		//         "MHDReader::readLevel | Error: Incorrect number of components. %i != %i",
		//         C, numComps);

		// READ DOMAIN SIZES
		int domainSize[2];
		auto att_domain = H5Aopen(*a_file, "domain", H5P_DEFAULT);
		assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0);
		assert(H5Aclose(att_domain) >= 0);

		int rCoord = CUBED_SPHERE_SHELL_RADIAL_COORD;
		int thetaCoord = (rCoord + 1) % 3;
		int phiCoord = (rCoord + 2) % 3;

		// these values are backwards in the data file
		Point sizes = Point::Ones();
		sizes[thetaCoord] = domainSize[0];
		sizes[phiCoord] = domainSize[1];

		return Box(sizes);
	}

	template <typename... Args>
	inline Box readDomain(
			std::string a_filename,
			Args... a_params)
	{
		hid_t file;
		readFile(&file, a_filename, a_params...);
		Box domainBox = _readDomain(&file);
		barrier();
		assert(H5Fclose(file) >= 0);
		return domainBox;
	}

	template <typename T, unsigned int C, MemType MEM, typename... Args>
	inline void readDataset(hid_t *a_file,
													BoxData<T, C, MEM> &a_data,
													int a_index)
	{
		hid_t H5T_T;
		HDF5Handler::getH5DataType<T>(&H5T_T);
		auto ds_data = H5Dopen2(*a_file, ("/data" + std::to_string(a_index)).c_str(), H5P_DEFAULT);
		if (MEM == HOST)
		{
			assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.data()) >= 0);
		}
		else
		{
			int dataSize = a_data.size();
			T *rawData = (T *)proto_malloc<HOST>(dataSize);
			assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawData) >= 0);
			proto_memcpy<HOST, MEM>(rawData, a_data.data(), dataSize);
			proto_free<HOST>(rawData);
		}
		assert(H5Dclose(ds_data) >= 0);
	}

	template <typename T, unsigned int C, MemType MEM, typename... Args>
	inline void readData(
			BoxData<T, C, MEM> &a_data,
			std::string a_filename,
			Args... a_params)
	{
		hid_t file;
		readFile(&file, a_filename, a_params...);
		Box domainBox = _readDomain(&file);
		a_data.define(domainBox);
		readDataset(&file, a_data, 0);
		// barrier();
		assert(H5Fclose(file) >= 0);
	}

	template <typename T, unsigned int C, MemType MEM, typename... Args>
	inline void readData(
			std::vector<BoxData<T, C, MEM>> &a_data,
			std::string a_filename,
			Args... a_params)
	{
		hid_t file;
		readFile(&file, a_filename, a_params...);
		Box domainBox = _readDomain(&file);

		int numDatasets;
		auto att_domain = H5Aopen(file, "num_datasets", H5P_DEFAULT);
		assert(H5Aread(att_domain, H5T_STD_I32LE, &numDatasets) >= 0);
		assert(H5Aclose(att_domain) >= 0);

		a_data.clear();
		a_data.resize(numDatasets);
		for (int ii = 0; ii < numDatasets; ii++)
		{
			a_data[ii].define(domainBox);
			readDataset(&file, a_data[ii], ii);
		}
		// barrier();
		assert(H5Fclose(file) >= 0);
	}

	template <typename T, typename... Args>
	inline void readGeom(
			std::vector<T> &a_data,
			std::string a_filename,
			Args... a_params)
	{
		PROTO_ASSERT(DIM == 3,
								 "MHDReader::readGeom | Error: Invalid for DIM != 3");

		PR_TIME("MHDReader::readGeom");
		hid_t file;
		readFile(&file, a_filename, a_params...);

		hid_t H5T_T;
		HDF5Handler::getH5DataType<T>(&H5T_T);

		// READ DOMAIN SIZES
		int domainSize[2];
		auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
		assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0);
		assert(H5Aclose(att_domain) >= 0);

		// READ DATA
		int dataSize = domainSize[1];
		auto ds_data = H5Dopen2(file, "/geometry/dtheta", H5P_DEFAULT);

		// DEFINE VECTOR
		a_data.resize(dataSize);

		// READ DATA
		assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_data[0]) >= 0);

		// CLEAN UP
		assert(H5Dclose(ds_data) >= 0);
		barrier();
		assert(H5Fclose(file) >= 0);
	}

	template <typename T, typename... Args>
	inline void readTheta(
			std::vector<T> &a_data,
			std::string a_filename,
			Args... a_params)
	{
		PROTO_ASSERT(DIM == 3,
								 "MHDReader::readGeom | Error: Invalid for DIM != 3");

		PR_TIME("MHDReader::readGeom");
		hid_t file;
		readFile(&file, a_filename, a_params...);

		hid_t H5T_T;
		HDF5Handler::getH5DataType<T>(&H5T_T);

		// READ DOMAIN SIZES
		int domainSize[2];
		auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
		assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0);
		assert(H5Aclose(att_domain) >= 0);

		// READ DATA
		int dataSize = domainSize[1];
		auto ds_data = H5Dopen2(file, "/geometry/theta", H5P_DEFAULT);

		// DEFINE VECTOR
		a_data.resize(dataSize);

		// READ DATA
		assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_data[0]) >= 0);

		// CLEAN UP
		assert(H5Dclose(ds_data) >= 0);
		barrier();
		assert(H5Fclose(file) >= 0);
	}

	template <typename T, typename... Args>
	inline void readPhi(
			std::vector<T> &a_data,
			std::string a_filename,
			Args... a_params)
	{
		PROTO_ASSERT(DIM == 3,
								 "MHDReader::readGeom | Error: Invalid for DIM != 3");

		PR_TIME("MHDReader::readGeom");
		hid_t file;
		readFile(&file, a_filename, a_params...);

		hid_t H5T_T;
		HDF5Handler::getH5DataType<T>(&H5T_T);

		// READ DOMAIN SIZES
		int domainSize[2];
		auto att_domain = H5Aopen(file, "domain", H5P_DEFAULT);
		assert(H5Aread(att_domain, H5T_STD_I32LE, domainSize) >= 0);
		assert(H5Aclose(att_domain) >= 0);

		// READ DATA
		int dataSize = domainSize[0];
		auto ds_data = H5Dopen2(file, "/geometry/phi", H5P_DEFAULT);

		// DEFINE VECTOR
		a_data.resize(dataSize);

		// READ DATA
		assert(H5Dread(ds_data, H5T_T, H5S_ALL, H5S_ALL, H5P_DEFAULT, &a_data[0]) >= 0);

		// CLEAN UP
		assert(H5Dclose(ds_data) >= 0);
		barrier();
		assert(H5Fclose(file) >= 0);
	}
};
