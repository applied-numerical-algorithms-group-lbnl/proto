// DO NOT INCLUDE THIS FILE. IT IS FOR DOCUMENTATION PURPOSES ONLY
#pragma once

#include "Proto.H"

using namespace Proto;

// Define a device-ready function for computing the low level map X0->X,J
// This function can have any signature so long as the inputs and outputs
// are designed to function on an accelerator. Point, Box, Array, and built
// in types like int and double will all work here.
namespace {
    ACCEL_DECORATION
    void f_map(Array<double, DIM>& a_X, double& J, const Array<double, DIM>& a_X0, int a_block)
    {
        // compute X, J from mapped coordinates X0 and the block 
    }
}

// This is the actual map operator
template<MemType MEM>
class MBMap_Template : public MBMapOp<MEM>
{
    public:

    // inherit constructors
    using MBMapOp<MEM>::MBMapOp;

    // this function must be defined. The arguments are both NODE CENTERED
    // and should have the same box.
    inline void apply(
            BoxData<double, DIM, MEM>& a_X,
            BoxData<double, 1, MEM>& a_J)
    {
        //inside MBMapOp the following context and API is available:
        auto& dx = this->dx(); // grid spacing (mapped)
        auto X0 = this->X(a_X.box(), dx); // get the mapped coordinates
        auto block = this->block(); // get the current mapped block

        // this is a sample implementation. X and J can be generated any way you like.
        forallInPlace(
        [] PROTO_LAMBDA (Var<double, DIM, MEM>& a_X, Var<double>& a_J, Var<double, DIM, MEM>& a_X0, int a_block)
        {
            Array<double, DIM> X0;
            for (int ii = 0; ii < DIM; ii++) { X0[ii] = a_X0(ii); }
            Array<double, DIM> X;
            double J;
            f_map(X, J, X0, a_block); // the function defined at the top of this file
            for (int ii = 0; ii < DIM; ii++) { a_X(ii) = X[ii]; }
            a_J(0) = J;
        }, a_X, a_J, X0, block
        );
    }
   
    // Computes the map for a single array. This function must also be defined to allow for analytic mapping. 
    // The following is a sample implementation which reuses the f_map function
    inline Array<double, DIM> apply(const Array<double, DIM>& a_X)
    {
        auto block = this->block();
        Array<double, DIM> X;
        double J; //unused
        f_map(X, J, a_X, block);
        return X;
    }
};

