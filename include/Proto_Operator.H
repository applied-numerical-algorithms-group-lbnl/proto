#pragma once

#include "Proto_BoxData.H"
//#include "AMRState.H"

namespace Operator {

    using namespace Proto;    

    /*
    template<unsigned int C>
    inline void faceDivergence(PatchState<C, 0>& a_div,
            const PatchState<0,C>& a_flux,
            double a_dx);
    */
    /// 4th Order Convolution
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> convolve(
            const BoxData<T,C,D,E>& a_ctr,
            const BoxData<T,C,D,E>& a_2nd);
    
    /// 4th Order Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void convolve(
            BoxData<T,C,D,E>& a_avg,
            const BoxData<T,C,D,E>& a_ctr,
            const BoxData<T,C,D,E>& a_2nd);

    /// 4th Order Deconvolution
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> deconvolve(
            const BoxData<T,C,D,E>& a_avg,
            const BoxData<T,C,D,E>& a_2nd);
    
    /// 4th Order Deconvolution (In Place)
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void deconvolve(
            BoxData<T,C,D,E>& a_ctr,
            const BoxData<T,C,D,E>& a_avg,
            const BoxData<T,C,D,E>& a_2nd);

    /// 4th Order Face Convolution
    /**
        Compute the 4th order accurate average of a quantity over a face. 
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> convolveFace(
            const BoxData<T,C,D,E>& a_ctr,
            const BoxData<T,C,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a face. 
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void convolveFace(
            BoxData<T,C,D,E>& a_avg,
            const BoxData<T,C,D,E>& a_ctr,
            const BoxData<T,C,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Deconvolution
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> deconvolveFace(
            const BoxData<T,C,D,E>& a_avg,
            const BoxData<T,C,D,E>& a_2nd,
            int a_dir);
   
    /// 4th Order Face Deconvolution (In Place)
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void deconvolveFace(
            BoxData<T,C,D,E>& a_ctr,
            const BoxData<T,C,D,E>& a_avg,
            const BoxData<T,C,D,E>& a_2nd,
            int a_dir);
   
    /// 4th Order Cell Centered Product Convolution (In Place)
    /**
        Compute the 4th order accurate value of: <p*q> == P(<p>, <q>, p_2, q_2)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> cellProduct(
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd);
    
    /// 4th Order Cell Centered Product Convolution (In Place)
    /**
        Compute the 4th order accurate value of: <p*q> == P(<p>, <q>, p_2, q_2)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void cellProduct(
            BoxData<T,C,D,E>& a_product,
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd);
    
    /// 4th Order Face Centered Product Convolution
    /**
        Compute the 4th order accurate value of: <p*q>_d == P(<p>d, <q>d, p_2d, q_2d)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> faceProduct(
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Centered Product Convolution (In Place)
    /**
        Compute the 4th order accurate value of: <p*q>_d == P(<p>d, <q>d, p_2d, q_2d)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void faceProduct(BoxData<T,C,D,E>& a_dest,
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Cell Centered Quotient Convolution
    /**
        Compute the 4th order accurate value of: <p/q> == Q(<p>, <q>, p_2, q_2)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> cellQuotient(
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd);
    
    /// 4th Order Cell Centered Quotient Convolution (In Place)
    /**
        Compute the 4th order accurate value of: <p/q> == Q(<p>, <q>, p_2, q_2)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void cellQuotient(
            BoxData<T,C,D,E>& a_quotient,
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd);
    
    /// 4th Order Face Centered Quotient Convolution
    /**
        Compute the 4th order accurate value of: <p/q>_d == Q(<p>d, <q>d, p_2d, q_2d)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> faceQuotient(
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd,
            int a_dir);

    /// 4th Order Face Centered Quotient Convolution (In Place)
    /**
        Compute the 4th order accurate value of: <p/q>_d == Q(<p>d, <q>d, p_2d, q_2d)
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void faceQuotient(
            BoxData<T,C,D,E>& a_quotient,
            const BoxData<T,C,D,E>& a_pavg,
            const BoxData<T,C,D,E>& a_qavg,
            const BoxData<T,C,D,E>& a_p2nd,
            const BoxData<T,C,D,E>& a_q2nd,
            int a_dir);

    /// 4th Order Cell to Face Interpolation
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> cellToFace(
            const BoxData<T,C,D,E>& a_cell,
            int a_dir);
    
    /// 4th Order Cell to Face Interpolation (In Place)
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void cellToFace(
            BoxData<T,C,D,E>& a_face,
            const BoxData<T,C,D,E>& a_cell,
            int a_dir);
    
    /// 4th Order Cell to Face Derivative (Undivided)
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline BoxData<T,C,D,E> diffCellToFace(
            const BoxData<T,C,D,E>& a_data,
            int a_dir);
    
    /// 4th Order Cell to Face Derivative (Undivided) (In Place)
    template<class T, unsigned int C, unsigned char D, unsigned char E>
    inline void diffCellToFace(
            BoxData<T,C,D,E>& a_deriv,
            const BoxData<T,C,D,E>& a_data,
            int a_dir);
}
#include "implem/Proto_OperatorI.H"
