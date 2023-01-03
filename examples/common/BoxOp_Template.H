#pragma once
#include "Proto.H"

#define NUMCOMPS 1
#define NUMAUX 0

using namespace Proto;

template<typename T, MemType MEM = MEMTYPE_DEFAULT>
class BoxOp_Template : public BoxOp<T, NUMCOMPS, NUMAUX, MEM>
{
    public:

    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(2);}
    
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Zeros();}
`
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() { return 2; }

    inline BoxOp_Template() : BoxOp<T, 1, 0, MEM>() {}
    inline BoxOp_Template(double a_dx) : BoxOp<T, 1, 0, MEM>(a_dx) {}
    inline BoxOp_Template(Array<double, DIM> a_dx) : BoxOp<T, 1, 0, MEM>(a_dx)
    {
        // Fill this out as needed if the operator has state other than that defined in BoxOp.H
        // or provided by auxiliary data
    }

    // The approximate spectral radius. Often needed if the operator is to be inverted using 
    // an iterative solver
    inline double spectralRadius() const {return 1;}
   
    //=============================================================================================
    // OPERATOR DEFINITIONS (BASIC) |
    //------------------------------+

    // Use only the Basic or Auxiliary functions to define your operator. The unused functions
    // can be deleted. Note that the grid spacing is available through this->dx().

    // Flux definition: Delete this if not needed.
    /*
    inline void flux(
            BoxData<T, NUMCOMPS>& a_flux,
            const BoxData<T, NUMCOMPS>& a_state,
            int a_dir) const
    {
    }
    */
    
    // Source definition: Delete this if not needed.
    /*
    inline void source(
            BoxData<T, NUMCOMPS>& a_source,
            const BoxData<T, NUMCOMPS>& a_state) const
    {
    }
    */
    
    // Diag definition: Delete this if not needed.
    /*
    inline void diag(
            BoxData<T, NUMCOMPS>& a_diag,
            const BoxData<T, NUMCOMPS>& a_state) const
    {
    }
    */
    
    // Apply definition: Delete this if not needed.
    // Only use this if you need additional control over the way the operator is computed
    // otherwise, this function will be inferred from the flux, source, and diag functions
    /*
    inline void operator()(
            BoxData<T, NUMCOMPS>&                   a_output,
            Array<BoxData<T, NUMCOMPS>, DIM>&  a_fluxes,
            const BoxData<T, NUCOMPS>&              a_state,
            T                                       a_scale = 1.0) const
    {
    }
    */
    
    //=============================================================================================
    // OPERATOR DEFINITIONS (AUX) |
    //----------------------------+

    // THESE ARE EXPERIMENTAL AT BEST. If you're really sure you need these, please let me know
    // so I can support your effort - CLG (cgebhart@lbl.gov)

    // Use only the Basic or Auxiliary functions to define your operator. The unused functions
    // can be deleted. Note that the grid spacing is available through this->dx().
    
    // Flux definition: Delete this if not needed.
    /*
    inline void flux(
            BoxData<T, NUMCOMPS>& a_flux,
            const BoxData<T, NUMCOMPS>& a_state,
            const BoxData<T, NUMAUX>& a_aux,
            int a_dir) const
    {
    }
    */
    
    // Source definition: Delete this if not needed.
    /*
    inline void source(
            BoxData<T, NUMCOMPS>& a_source,
            const BoxData<T, NUMCOMPS>& a_state,
            const BoxData<T, NUMAUX>& a_aux) const
    {
    }
    */
    
    // Diag definition: Delete this if not needed.
    /*
    inline void diag(
            BoxData<T, NUMCOMPS>& a_diag,
            const BoxData<T, NUMCOMPS>& a_state,
            const BoxData<T, NUMAUX>& a_aux) const
    {
    }
    */
    
    // Apply definition: Delete this if not needed.
    // Only use this if you need additional control over the way the operator is computed
    // otherwise, this function will be inferred from the flux, source, and diag functions
    /*
    inline void operator()(
            BoxData<T, NUMCOMPS>&                   a_output,
            Array<BoxData<T, NUMCOMPS>, DIM>&  a_fluxes,
            const BoxData<T, NUCOMPS>&              a_state,
            const AuxData&                          a_aux,
            T                                       a_scale = 1.0) const
    {
    }
    */
};
