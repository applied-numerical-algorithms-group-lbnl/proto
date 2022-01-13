#pragma once
#ifndef PROTO_HDF5_H
#define PROTO_HDF5_H
#include "Proto_LevelBoxData.H"
#include "Proto_AMRData.H"
#include "hdf5.h"

namespace Proto
{
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
        template<unsigned int C>
            void writeHeader(
                    hid_t* a_file,
                    int a_numLevels,
                    std::vector<std::string> a_varNames);
        /// Write Level Box Data
        /**
            Write's a LevelBoxData to an HDF5 file which is compatible with Chombo's HDF5 I/O
            as well as VisIt.
            
            The filename and following variadic arguments are used in the same way as 
            printf-style functions. For example:

            writeLevel(data, {"rho", "vx", "vy"}, dx, "Euler_T%i.hdf5", 0) -> writes to "Euler_T0.hdf5"
            
            Compatible with MPI if PR_MPI=TRUE.
        */
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeLevel(
                        std::vector<std::string> a_varnames,
                        array<double, DIM> a_dx,
                        const LevelBoxData<T, C, MEM, CTR>& a_data,
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
                        array<double, DIM> a_dx,
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
      
        /// Write AMR Hierarchy 
        template<typename T, unsigned int C, MemType MEM, Centering CTR, typename... Args>
        inline void writeAMRData(
                        std::vector<std::string> a_varnames,
                        array<double, DIM> a_dx,
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
                        array<double, DIM> a_dx,
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
        
        /// Write Box Data
        /**
            Writes a single BoxData to an HDF5 file compatible with Chombo's HDF5 I/O and VisIt.
            
            The filename and following variadic arguments are used in the same way as 
            printf-style functions. For example:

            writePatch(data, {"rho", "vx", "vy"}, dx, "Euler_T%i.hdf5", 0) -> writes to "Euler_T0.hdf5"
            
            Written patches are shifted such that their domain is at (0,...,0) by necessity and,
            if they were originally part of a LevelBoxData, have no notion of ghost cells or
            their previous location in the Box layout.

            Ghost cells ARE printed, but appear with the rest of the data (they are not hidden by
            default in VisIt).
            
            This function should only be used for simple debugging needs wherein it is inconvenient to
            print an entire LevelBoxData.
        
            TODO: Does not currently work with MPI.
        */
        template<typename T, unsigned int C, MemType MEM, typename... Args>
        inline void writePatch(
                        std::vector<std::string> a_varnames,
                        array<double, DIM> a_dx,
                        const BoxData<T, C, MEM>& a_data,
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
                        array<double, DIM> a_dx,
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
        private:
      
        hid_t H5T_PROTO_POINT;
        hid_t H5T_PROTO_BOX; 
        hid_t H5T_PROTO_STRING;

        template<typename T, unsigned int C, MemType MEM, Centering CTR>
        inline void addLevel(hid_t* a_file,
            const LevelBoxData<T, C, MEM, CTR>& a_data,
            array<double, DIM>& a_dx0,
            Point a_refRatio,
            int a_level);
        
        template<typename T>
        inline void getH5DataType(hid_t* a_type) {}
        
    }; // end class HDF5Handler
#include "implem/Proto_HDF5Implem.H"
} // end namespace Proto
#endif //end of include guard
