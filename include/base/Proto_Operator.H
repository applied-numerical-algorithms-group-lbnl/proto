#pragma once
#ifndef _PROTO_OPERATOR_H_
#define _PROTO_OPERATOR_H_
#include "Proto_HDF5.H"
#include "Proto.H"

#include "Proto_BoxData.H"
#include "Proto_HDF5.H"
#include "Proto_LevelBoxData.H"

namespace Operator {

    using namespace Proto;    

    /// 4th Order Convolution
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Distinguishes between 4th and 2nd order inputs.
         
        \param a_ctr    4th  order cell centered quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _convolve(
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd);
    
    /// 4th Order Convolution
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Equal to _convolve(ctr, ctr).

        \param a_ctr    4th  order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> convolve(
            const BoxData<T,C,MEM,D,E>& a_ctr);
    
    /// 4th Order Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Distinguishes between 4th and 2nd order inputs.
       
        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _convolve(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd);
    
    /// 4th Order Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a cell. 
        Equal to _convolve(avg, ctr, ctr).

        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void convolve(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr);
    
    /// 4th Order Convolution (Level, In Place)
    /**
        Compute the 4th order accurate average of a quantity over a level. 
        Distinguishes between 4th and 2nd order inputs.
        TODO: implement appropriate behavior for non-cell centered LevelData
         
        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void _convolve(
            LevelBoxData<T,C,MEM,CTR>& a_avg,
            const LevelBoxData<T,C,MEM,CTR>& a_ctr,
            const LevelBoxData<T,C,MEM,CTR>& a_2nd);
    /// 4th Order Convolution (Level, In Place)
    /**
        Compute the 4th order accurate average of a quantity over a level. 
        Equal to _convolve(avg, ctr, ctr).
        TODO: implement appropriate behavior for non-cell centered LevelData
        
        \param a_avg    4th  order cell averaged quantity (output)
        \param a_ctr    4th  order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void convolve(
            LevelBoxData<T,C,MEM,CTR>& a_avg,
            const LevelBoxData<T,C,MEM,CTR>& a_ctr);
    
    /// Initialize and Convolve
    template<class T, unsigned int C, MemType MEM, Centering CTR,
        typename Func, typename... Srcs>
    inline void initConvolve(LevelBoxData<T,C,MEM,CTR>& a_data, Func& a_func, Srcs... a_srcs);
    
#ifdef PR_AMR
    /// Initialize and Convolve (AMR)
    template<class T, unsigned int C, MemType MEM, Centering CTR,
        typename Func, typename... Srcs>
    inline void initConvolve(AMRData<T,C,MEM,CTR>& a_data, double a_dx0, Func& a_func, Srcs... a_srcs);
#endif

    /// 4th Order Deconvolution
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Distinguishes between 4th and 2nd order inputs

        \param a_avg    4th  order cell averaged quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _deconvolve(
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd);
    
    /// 4th Order Deconvolution
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Equal to _deconvolve(avg, avg).
            
        \param a_avg    4th  order cell averaged quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> deconvolve(
            const BoxData<T,C,MEM,D,E>& a_avg);
    
    /// 4th Order Deconvolution (In Place)
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Distinguishes between 4th and 2nd order inputs
        
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _deconvolve(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd);

    /// 4th Order Deconvolution (In Place)
    /**
        Compute the 4th order accurate cell centered value of a quantity from its cell average.
        Equal to _deconvolve(ctr, avg, avg)
        
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void deconvolve(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg);
    
    /// 4th Order Deconvolution (Level, In Place)
    /**
        Compute the 4th order accurate cell centered quantity over a level. 
        Distinguishes between 2nd and 4th order inputs.
        TODO: implement appropriate behavior for non-cell centered LevelData
         
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
        \param a_2nd    2nd+ order cell centered quantity
    */
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void _deconvolve(
            LevelBoxData<T,C,MEM,CTR>& a_ctr,
            const LevelBoxData<T,C,MEM,CTR>& a_avg,
            const LevelBoxData<T,C,MEM,CTR>& a_2nd);
    /// 4th Order Deconvolution (Level, In Place)
    /**
        Compute the 4th order accurate cell centered quantity over a level. 
        Equal to _deconvolve(ctr, avg, avg);
        TODO: implement appropriate behavior for non-cell centered LevelData
         
        \param a_ctr    4th  order cell centered quantity (output)
        \param a_avg    4th  order cell averaged quantity
    */
    template<class T, unsigned int C, MemType MEM, Centering CTR>
    inline void deconvolve(
            LevelBoxData<T,C,MEM,CTR>& a_ctr,
            const LevelBoxData<T,C,MEM,CTR>& a_avg);
    /// 4th Order Face Convolution
    /**
        Compute the 4th order accurate average of a quantity over a face.
        Distinguishes between 2nd and 4th order inputs.

        \param a_ctr    4th  order face centered quantity
        \param a_2nd    2nd+ order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _convolveFace(
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Convolution
    /**
        Compute the 4th order accurate average of a quantity over a face.
        Equal to _convolveFace(ctr, ctr, dir);

        \param a_ctr    4th  order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> convolveFace(
            const BoxData<T,C,MEM,D,E>& a_ctr,
            int a_dir);
    
    /// 4th Order Face Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a face. 
        Distinguishes between 2nd and 4th order inputs.

        \param a_avg    4th  order face averaged quantity (output)
        \param a_ctr    4th  order face centered quantity
        \param a_2nd    2nd+ order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _convolveFace(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Convolution (In Place)
    /**
        Compute the 4th order accurate average of a quantity over a face. 
        Equal to _convolveFace(avg, ctr, ctr, dir)
        
        \param a_avg    4th  order face averaged quantity (output)
        \param a_ctr    4th  order face centered quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void convolveFace(
            BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_ctr,
            int a_dir);
    
    /// 4th Order Face Deconvolution
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Distinguishes between 2nd and 4th order inputs.

        \param a_avg    4th  order face averaged quantity
        \param a_2nd    2nd+ order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _deconvolveFace(
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Deconvolution
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Equal to _deconvolveFace(avg, avg, dir)
        
        \param a_avg    4th  order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> deconvolveFace(
            const BoxData<T,C,MEM,D,E>& a_avg,
            int a_dir);
   
    /// 4th Order Face Deconvolution (In Place)
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Distinguishes between 2nd and 4th order inputs.

        \param a_ctr    4th  order face centered quantity (output)
        \param a_avg    4th  order face averaged quantity
        \param a_2nd    2nd+ order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _deconvolveFace(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg,
            const BoxData<T,C,MEM,D,E>& a_2nd,
            int a_dir);
    
    /// 4th Order Face Deconvolution (In Place)
    /**
        Compute the 4th order accurate face centered value of a quantity from its face average.
        Equal to _deconvolveFace(ctr, avg, avg, dir)
         
        \param a_ctr    4th  order face centered quantity (output)
        \param a_avg    4th  order face averaged quantity
        \param a_dir    Face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void deconvolveFace(
            BoxData<T,C,MEM,D,E>& a_ctr,
            const BoxData<T,C,MEM,D,E>& a_avg,
            int a_dir);
   
    /// 4th Order Cell Averaged Product
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Distinguishes between 2nd and 4th order inputs.
    
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _cellProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Product
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Equal to _cellProduct(pavg, qavg, pavg, qavg)
    
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> cellProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Cell Averaged Product (In Place)
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Distinguishes between 2nd and 4th order inputs.
        
        \param a_prod   4th  order cell averaged product P*Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _cellProduct(
            BoxData<T,C,MEM,D,E>& a_prod,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Product (In Place)
    /**
        Compute the 4th order accurate cell averaged product rule <P>*<Q> 
        Equal to _cellProduct(prod, pavg, qavg, pavg, qavg)

        \param a_prod   4th  order cell averaged product P*Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void cellProduct(
            BoxData<T,C,MEM,D,E>& a_product,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Face Averaged Product
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Distinguishes between 2nd and 4th order inputs.
        
        \param a_prod   4th  order face averaged product P*Q (output)
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal 
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _faceProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Product
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Equal to _faceProduct(pavg, qavg, pavg, qavg, dir)
        
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> faceProduct(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);
    
    /// 4th Order Face Averaged Product (In Place)
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Distinguishes between 2nd and 4th order inputs.

        \param a_prod   4th  order face averaged product P*Q (output)
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _faceProduct(
            BoxData<T,C,MEM,D,E>& a_prod,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Product (In Place)
    /**
        Compute the 4th order accurate face averaged product rule <P>_d * <Q>_d 
        Equal to _faceProduct(prod, pavg, qavg, pavg, qavg, dir)
        
        \param a_prod   4th  order face averaged product P*Q (output)
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void faceProduct(
            BoxData<T,C,MEM,D,E>&       a_prod,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);
    
    /// 4th Order Cell Averaged Quotient
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Distinguishes between 2nd and 4th order inputs.

        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _cellQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Quotient
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Equal to _cellQuotient(pavg, qavg, pavg, qavg)

        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> cellQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Cell Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Distinguishes between 2nd and 4th order inputs.

        \param a_quot   4th  order cell averaged quotient P / Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
        \param a_p2nd   2nd+ order cell centered quantity P
        \param a_q2nd   2nd+ order cell centered quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _cellQuotient(
            BoxData<T,C,MEM,D,E>&       a_quot,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd);
    
    /// 4th Order Cell Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate cell averaged quotient rule <P> / <Q> 
        Equal to _cellQuotient(qout, pavg, qavg, pavg, qavg)

        \param a_quot   4th  order cell averaged quotient P / Q (output)
        \param a_pavg   4th  order cell averaged quantity P
        \param a_qavg   4th  order cell averaged quantity Q
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void cellQuotient(
            BoxData<T,C,MEM,D,E>&       a_quot,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg);
    
    /// 4th Order Face Averaged Quotient
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Distinguishes between 2nd and 4th order inputs.

        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> _faceQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir);
    
    /// 4th Order Face Averaged Quotient
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Equal to _faceQuotient(pavg, qavg, pavg, qavg, dir)
            
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline BoxData<T,C,MEM,D,E> faceQuotient(
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);

    /// 4th Order Face Averaged Quotient (In Place)
    /**
        Compute the 4th order accurate face averaged quotient rule <P>_d / <Q>_d 
        Distinguishes between 2nd and 4th order inputs.
                
        \param a_quot   4th  order face averaged quotient P / Q
        \param a_pavg   4th  order face averaged quantity P
        \param a_qavg   4th  order face averaged quantity Q
        \param a_p2nd   2nd+ order face centered quantity P
        \param a_q2nd   2nd+ order face centered quantity Q
        \param a_dir    face normal
    */
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void _faceQuotient(
            BoxData<T,C,MEM,D,E>& a_quotient,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            const BoxData<T,C,MEM,D,E>& a_p2nd,
            const BoxData<T,C,MEM,D,E>& a_q2nd,
            int a_dir); 
    
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    inline void faceQuotient(
            BoxData<T,C,MEM,D,E>& a_quotient,
            const BoxData<T,C,MEM,D,E>& a_pavg,
            const BoxData<T,C,MEM,D,E>& a_qavg,
            int a_dir);
    
    /// 4th Order Face Averaged Tensor Quotient 
    /**
         Returns the 4th order accurate face averaged quotient rule <(1/QR)*QL>
         where QL is a tensor,and QR is a scalar. 
         Distinguishes between 2nd and 4th order inputs.
        
        \param a_ql4    4th  order face averaged tensor quantity QL
        \param a_qr4    4th  order face averaged scalar quantity QR
        \param a_ql2    2nd+ order face centered quantity QL
        \param a_qr2    2nd+ order face centered quantity QR
        \param a_dir    face normal
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
    BoxData<T,C,MEM,D,E>
    _faceTensorQuotient(
            const BoxData<T,C,MEM,D,E>& a_ql4,
            const BoxData<T,1,MEM,1,1>& a_qr4,
            const BoxData<T,C,MEM,D,E>& a_ql2,
            const BoxData<T,1,MEM,1,1>& a_qr2,
            int a_dir);
 /// 4th Order Face Averaged Tensor Product 
    /**
         Returns the 4th order accurate face averaged quotient rule <QR*QL>
         where QL is a tensor,and QR is a scalar. 
         Distinguishes between 2nd and 4th order inputs.
        
        \param a_ql4    4th  order face averaged tensor quantity QL
        \param a_qr4    4th  order face averaged scalar quantity QR
        \param a_ql2    2nd+ order face centered quantity QL
        \param a_qr2    2nd+ order face centered quantity QR
        \param a_dir    face normal
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
    BoxData<T,C,MEM,D,E>
    _faceTensorProduct(
            const BoxData<T,C,MEM,D,E>& a_ql4,
            const BoxData<T,1,MEM,1,1>& a_qr4,
            const BoxData<T,C,MEM,D,E>& a_ql2,
            const BoxData<T,1,MEM,1,1>& a_qr2,
            int a_dir);
    /// 4th Order Cell Averaged Tensor Quotient 
    /**
         Returns the 4th order accurate face averaged quotient rule <(1/QR) QL>
         where QL is a tensor,and QR is a scalar. 
         Distinguishes between 2nd and 4th order inputs.
        
        \param a_ql4    4th  order cell averaged tensor quantity QL
        \param a_qr4    4th  order cell averaged scalar quantity QR
        \param a_ql2    2nd+ order cell centered quantity QL
        \param a_qr2    2nd+ order cell centered quantity QR
    */
    template<class T, unsigned int C, unsigned char D, unsigned char E, MemType MEM>
    BoxData<T,C,MEM,D,E>
    _cellTensorQuotient(                    
            const BoxData<T,C,MEM,D,E>& a_ql4,
            const BoxData<T,1,MEM,1,1>& a_qr4,
            const BoxData<T,C,MEM,D,E>& a_ql2,
            const BoxData<T,1,MEM,1,1>& a_qr2);

    /// Cofactor Vector On Face
    template <class T, MemType MEM>
    BoxData<T,DIM,MEM>
    cofactor(
            const BoxData<T,DIM,MEM>& a_X,
            int a_dir);

    /// Compute Matrix of Cofactors
    template<class T, MemType MEM>
    inline BoxData<T,DIM,MEM,DIM>
    _cofactorMatrix(
            const Array<BoxData<T,DIM,MEM>, DIM>& a_NTFaces,
            int a_dirNorm);

    template <class T, MemType MEM>
    inline BoxData<T,1,MEM>
    jacobian(
            const BoxData<T,DIM,MEM>& a_X,
            const Array<BoxData<T,DIM,MEM>, DIM>& a_NT);

    template<class T, MemType MEM>
    BoxData<T,3,MEM>
    _edgeCrossProduct3D(
            const BoxData<T,3,MEM>& a_ql4,
            const BoxData<T,3,MEM>& a_qr4,
            const BoxData<T,3,MEM>& a_ql2,
            const BoxData<T,3,MEM>& a_qr2,
            int a_edgedir);

    template<class T, unsigned int CL,unsigned char DL,
           unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
    BoxData<T,CL,MEM,DR,E>
    _faceMatrixProductAB(
            const BoxData<T,CL,MEM,DL,E>& a_ql4,
            const BoxData<T,CR,MEM,DR,E>& a_qr4,
            const BoxData<T,CL,MEM,DL,E>& a_ql2,
            const BoxData<T,CR,MEM,DR,E>& a_qr2,
            int a_dir);

    template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
    BoxData<T,DL,MEM,DR,E>
    _faceMatrixProductATB(
            const BoxData<T,CL,MEM,DL,E>& a_ql4,
            const BoxData<T,CR,MEM,DR,E>& a_qr4,
            const BoxData<T,CL,MEM,DL,E>& a_ql2,
            const BoxData<T,CR,MEM,DR,E>& a_qr2,
            int a_dir);
    
    template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
    BoxData<T,CL,MEM,CR,E>
    _faceMatrixProductABT(
            const BoxData<T,CL,MEM,DL,E>& a_ql4,
            const BoxData<T,CR,MEM,DR,E>& a_qr4,
            const BoxData<T,CL,MEM,DL,E>& a_ql2,
            const BoxData<T,CR,MEM,DR,E>& a_qr2,
            int a_dir);
  
    template<class T, unsigned int CL,unsigned char DL,
           unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
    BoxData<T,CL,MEM,DR,E>
    _matrixProductAB2(
            const BoxData<T,CL,MEM,DL,E>& a_ql,
            const BoxData<T,CR,MEM,DR,E>& a_qr);

    template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
    BoxData<T,DL,MEM,DR,E>
    _matrixProductATB2(
            const BoxData<T,CL,MEM,DL,E>& a_ql,
            const BoxData<T,CR,MEM,DR,E>& a_qr);
    
    template<class T, unsigned int CL, unsigned char DL,
         unsigned int CR, unsigned char DR, unsigned char E, MemType MEM>
    BoxData<T,CL,MEM,CR,E>
    _matrixProductABT2(
            const BoxData<T,CL,MEM,DL,E>& a_ql,
            const BoxData<T,CR,MEM,DR,E>& a_qr);
  
    template<class T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
    BoxData<T,C,MEM,D,E>
    _product1D(
            const BoxData<T,C,MEM,D,E>& a_ql4,
            const BoxData<T,C,MEM,D,E>& a_qr4,
            const BoxData<T,C,MEM,D,E>& a_ql2,
            const BoxData<T,C,MEM,D,E>& a_qr2,
            int a_dir);
    
    template<class T,unsigned int C, MemType MEM>
    BoxData<T,C,MEM,DIM>
    _faceGradient(
                  const BoxData<T,C,MEM>& phi,
                  int a_dir,
                  int a_order=4);
  
  template<class T,unsigned int CFLUX,unsigned int CPRIM, MemType MEM>
  BoxData<T,CFLUX,MEM>
  MHDSphericalFlux(                     
                   const BoxData<T,CPRIM,MEM>& a_prim4,
                   const BoxData<T,CPRIM,MEM>& a_prim2,
                   const BoxData<T,DIM,MEM,DIM>& a_A4,
                   const BoxData<T,DIM,MEM,DIM>& a_A2,
                   const BoxData<T,1,MEM>& a_DrDetA4,                    
                   const BoxData<T,1,MEM>& a_DrDetA2,
                   const BoxData<T,DIM,MEM>& a_DrAdjA4,                    
                   const BoxData<T,DIM,MEM>& a_DrAdjA2,
                   const T& a_gamma,
                   int a_dir);
  template<class T,unsigned int C, MemType MEM>
  BoxData<T,C,MEM,DIM>
  _faceGradxPhi(                     
               const BoxData<T,C,MEM>& a_phi4,
               const BoxData<T,C,MEM>& a_phi2,
               const BoxData<T,DIM,MEM,DIM>& a_NFace4,
               const BoxData<T,DIM,MEM,DIM>& a_NFace2,
               const BoxData<T,1,MEM>& a_JFace4,
               const BoxData<T,1,MEM>& a_JFace2,
               int a_dir);
#include "implem/Proto_OperatorImplem.H"
}
#endif
