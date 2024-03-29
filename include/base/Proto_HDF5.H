#pragma once
#ifndef PROTO_HDF5_H
#define PROTO_HDF5_H
#include "Proto_BoxData.H"
#ifdef PR_HDF5
#include "hdf5.h"
#endif

namespace Proto
{
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class LevelBoxData;
    
#ifdef PR_AMR
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class AMRData;
    template<typename T, unsigned int C, MemType MEM>
    class LevelFluxRegister;
#endif
#ifdef PR_MMB
    template<MemType MEM>
    class MBLevelMap;
    template<typename T, unsigned int C, MemType MEM, Centering CTR>
    class MBLevelBoxData;
#endif

#ifdef PR_HDF5
    class HDF5Handler
    {
        public:

        /// Default Constructor
        inline HDF5Handler();
      
        /// Read Level Box Data
        /**
            Read a LevelBoxData from an appropriately constructed HDF5 file. This function 
            should be compatible with files created using Chombo's HDF5 I/O system.
            
            The input LevelBoxData is REDEFINED when calling this function, and a new
            DisjointBoxLayout is built. Because of this, knowledge of the way in which
            the input data was constructed is not needed. 

            The filename and following variadic arguments are used in the same way as 
            printf-style functions. For example:

            readLevel(data, "filename_%i.hdf5", 42) -> read from "filename_42.hdf5"
        
            Compatible with MPI if PR_MPI=TRUE.
        */ 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void readLevel(LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);

        /// Write Header
        template<unsigned int C, Centering CTR>
            void writeHeader(
                    hid_t* a_file,
                    int a_numLevels,
                    std::vector<std::string> a_varNames);
        
        /// Write Box Data
        /**
            Writes a single BoxData to an HDF5 file compatible with Chombo's HDF5 I/O and VisIt.
            
            The filename and following variadic arguments are used in the same way as 
            printf-style functions. For example:

            writePatch(data, {"rho", "vx", "vy"}, dx, "Euler_T%i.hdf5", 0) -> writes to "Euler_T0.hdf5"
            
            Ghost cells ARE printed, but appear with the rest of the data (they are not hidden by
            default in VisIt).
        */
        template<typename T, unsigned int C, typename... Args>
        inline void writePatch(
                        std::vector<std::string> a_varnames,
                        Array<double, DIM> a_dx,
                        const BoxData<T, C, HOST>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, typename... Args>
        inline void writePatch(
                        std::vector<std::string> a_varnames,
                        Array<double, DIM> a_dx,
                        const BoxData<T, C, DEVICE>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Box Data (Scalar DX)
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writePatch(
                        std::vector<std::string> a_varnames,
                        double a_dx,
                        const BoxData<T, C, MEM>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writePatch(
                        Array<double, DIM> a_dx,
                        const BoxData<T, C, MEM>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Box Data (Scalar DX)
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writePatch(
                        double a_dx,
                        const BoxData<T, C, MEM>& a_data,
                        std::string a_filename,
                        Args... a_params);
        /// Write Box Data (No DX)
        /**
            A convenient overload. Grid spacing defaults to 1.
        */ 
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writePatch(
                        std::vector<std::string> a_varnames,
                        const BoxData<T, C, MEM>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Box Data (No DX, Varnames)
        /**
            A convenient overload. Grid spacing defaults to 1 and variable names are generated
            automatically.
        */ 
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writePatch(const BoxData<T, C, MEM>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Level Box Data
        /**
            Write's a LevelBoxData to an HDF5 file which is compatible with Chombo's HDF5 I/O
            as well as VisIt.
            
            The filename and following variadic arguments are used in the same way as 
            printf-style functions. For example:

            writeLevel(data, {"rho", "vx", "vy"}, dx, "Euler_T%i.hdf5", 0) -> writes to "Euler_T0.hdf5"
            
            Compatible with MPI if PR_MPI=TRUE.
        */
        template<typename T, unsigned int C, Centering CTR, typename... Args>
        inline void writeLevel(
                        std::vector<std::string> a_varnames,
                        Array<double, DIM> a_dx,
                        const LevelBoxData<T, C, HOST, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        template<typename T, unsigned int C, Centering CTR, typename... Args>
        inline void writeLevel(
                        std::vector<std::string> a_varnames,
                        Array<double, DIM> a_dx,
                        const LevelBoxData<T, C, DEVICE, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Level Box Data (scalar DX)
        /**
            A convenient overload. Grid spacing is isotropic.
        */
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeLevel(
                        std::vector<std::string> a_varnames,
                        double a_dx,
                        const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Level Box Data (No Varnames)
        /**
            A convenient overload.
        */
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeLevel(
                        Array<double, DIM> a_dx,
                        const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Level Box Data (No Varnames, scalar DX)
        /**
            A convenient overload. Grid spacing is isotropic.
        */
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeLevel(
                        double a_dx,
                        const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Level Box Data (No DX)
        /**
            A convenient overload. Grid spacing defaults to 1.
        */ 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeLevel(
                        std::vector<std::string> a_varnames,
                        const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Level Box Data (No DX, Varnames)
        /**
            A convenient overload. Grid spacing defaults to 1 and variable names are generated
            automatically.
        */ 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeLevel(
                        const LevelBoxData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
#ifdef PR_AMR 
        /// Write AMR Hierarchy 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeAMRData(
                        std::vector<std::string> a_varnames,
                        Array<double, DIM> a_dx,
                        const AMRData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write AMR Hierarchy 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeAMRData(
                        std::vector<std::string> a_varnames,
                        double a_dx,
                        const AMRData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write AMR Hierarchy 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeAMRData(
                        Array<double, DIM> a_dx,
                        const AMRData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write AMR Hierarchy 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeAMRData(
                        double a_dx,
                        const AMRData<T, C, MEM, CTR>& a_data,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Flux Register
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writeFluxRegister(
                        std::vector<std::string> a_varnames,
                        Array<double, DIM> a_dx,
                        LevelFluxRegister<T, C, MEM>& a_register,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Flux Register
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writeFluxRegister(
                        Array<double, DIM> a_dx,
                        LevelFluxRegister<T, C, MEM>& a_register,
                        std::string a_filename,
                        Args... a_params);
        
        /// Write Flux Register
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writeFluxRegister(
                        LevelFluxRegister<T, C, MEM>& a_register,
                        std::string a_filename,
                        Args... a_params);
#endif
#ifdef PR_MMB
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeMBLevel(
                std::vector<std::string> a_varNames,
                const MBLevelMap<MEM>& a_map,
                const MBLevelBoxData<T, C, MEM, CTR>& a_data,
                std::string a_filename,
                Args... a_params);

        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeMBLevel(
                std::vector<std::string> a_varNames,
                const MBLevelBoxData<T, C, MEM, CTR>& a_data,
                std::string a_filename,
                Args... a_params);


        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeMBLevelBounds(
                std::vector<std::string> a_varNames,
                const MBLevelBoxData<T, C, MEM, CTR>& a_data,
                std::string a_filename,
                Args... a_params);

        template<typename T, unsigned int C, MemType MEM, Centering CTR>
        inline void addMBLevel(
                hid_t* a_file,
                const MBLevelBoxData<T, C, MEM, CTR>& a_data,
                Point a_refRatio,
                int a_level);
#endif
        inline void setTime(double a_time) { m_time = a_time; }
        inline void setTimestep(double a_dt) { m_dt = a_dt; }
        inline double& time() { return m_time; }
        inline double& dt() { return m_dt; }
        
        template<typename T>
        inline static void getH5DataType(hid_t* a_type) {}

        private:
      
      static inline hid_t H5T_PROTO_POINT();
      static inline hid_t H5T_PROTO_BOX(); 
      static inline hid_t H5T_PROTO_STRING();

        double m_time = 0.0;
        double m_dt = 1.0;

        template<typename T, unsigned int C, Centering CTR>
        inline void addLevel(hid_t* a_file,
            const LevelBoxData<T, C, DEVICE, CTR>& a_data,
            Array<double, DIM>& a_dx0,
            Point a_refRatio,
            int a_level);
        
        template<typename T, unsigned int C, Centering CTR>
        inline void addLevel(hid_t* a_file,
            const LevelBoxData<T, C, HOST, CTR>& a_data,
            Array<double, DIM>& a_dx0,
            Point a_refRatio,
            int a_level);
        
        template<typename T, unsigned int C, MemType MEM>
        inline void addPatch(hid_t* a_file,
            const BoxData<T, C, MEM>& a_data,
            Array<double, DIM>& a_dx0,
            Point a_ghost);
        
    }; // end class HDF5Handler
#include "implem/Proto_HDF5Implem.H"
#endif
} // end namespace Proto
#endif //end of include guard
