// DO NOT INCLUDE THIS FILE. IT IS FOR DOCUMENTATION PURPOSES ONLY
#pragma once
#include "Proto.H"

#define NSTATE
#define NAUX

using namespace Proto;

// The second template parameter is the MBMapOp template
// The third template parameter is the number of state variables
// The fourth template parameter is the number of auxiliary variables (usually 0)
template<typename T,
    template<MemType> typename MAP,
    MemType MEM = MEMTYPE_DEFAULT>
class MBBoxOp_Template : public MBBoxOp<T, MAP, NSTATE, NAUX, MEM>
{
    public:

    // inherit constructors
    using MBBoxOp<T, MAP, 1, 0, MEM>::MBBoxOp;
    
    // How many ghost cells does the operator need from the state variables
    inline static Point ghost() { return Point::Ones(1);}
    // How many ghost cells does the operator need from the auxiliary variables
    inline static Point auxGhost() { return Point::Ones(1);}
    
    // What is the intended order of accuracy of the operator
    inline static constexpr int order() const { return 4; }

    // An approximate value of the spectral radius is needed for iterative methods
    inline T spectralRadius() const {return 1.0; }

    // Init is called when the operator is constructed. This is a good place to cache
    // pieces of data like stencils and metric quantities. 
    void init(){}

    // If NAUX = 0, define this function
    inline void
    operator()(
        BoxData<T,NSTATE>& a_output,
        Array<BoxData<T,NSTATE>, DIM>& a_fluxes,
        const BoxData<T,NSTATE>& a_state,
        T                        a_scale = 1.0) const
    {
        // inside of a BoxOp, the following context is available
        Array<double, DIM> dx = this->dx(); //grid spacing (mapped coordinates)
        Box b = this->box(); // local domain (not including ghost cells)
        int rkStage = this->RKStage(); //stage of Runge-Kutta execution
        T t = this->time(); //timestep
        MBIndex ind = this->index(); //data index
        const MBDisjointBoxLayout& layout = this->layout(); //data layout
        const MBLevelMap<MAP,MEM>& map = this->map(); //map
        unsigned int block = this->block(); //mapped block
        const BoxData<double, DIM>& x = this->x(); //cached coordinates (also in map)
        const BoxData<double, 1>& J = this->jacobian(); //cached jacobian (also in map)
    }
    
    // If NAUX != 0, define this function
    inline virtual void
    operator()(
        BoxData<T,NSTATE>& a_output,
        Array<BoxData<T,NSTATE>, DIM>& a_fluxes,
        const BoxData<T,NSTATE>& a_state,
        const BoxData<T,NAUX>&   a_aux,
        T                        a_scale = 1.0) const {}
    
    // define as many additional helper functions as desired
    
    private:
    // BoxOps can have state if desired. This state can be initialized in the
    // init function which is called during construction of the BoxOp

};
