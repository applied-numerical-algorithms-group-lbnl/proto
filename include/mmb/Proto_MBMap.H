#pragma once
#ifndef _PROTO_MB_MAP_
#define _PROTO_MB_MAP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"
#include "Proto_Operator.H"

namespace Proto
{

#if 1

#define PROTO_MAP_START(name) \
    ACCEL_DECORATION \
    void proto_map_##name (double& x, double& y, double& z, \
            const double& X, const double& Y, const double& Z, unsigned int block) { \
        x = 0; y = 0; z = 0; \

#if DIM==3
#define PROTO_MAP_END(name) \
    } \
    struct struct_proto_map_##name { \
    template<MemType MEM = MEMTYPE_DEFAULT> \
    inline ACCEL_DECORATION \
    void operator() (Point& a_pt, Var<double, DIM, MEM>& a_x, \
            unsigned int a_block, \
            Point a_blockSize, \
            Centering a_centering) const \
    { \
        Array<double, DIM> dX; \
        for (int dir = 0; dir < DIM; dir++) { dX[dir] = 1.0/a_blockSize[dir]; } \
        Array<double, DIM> X0;\
        for (int dir = 0; dir < DIM; dir++) \
        { \
            if (a_centering == PR_NODE) \
            {\
                X0[dir] = a_pt[dir]*dX[dir]; \
            } else if (a_centering == PR_CELL) \
            {\
                X0[dir] = a_pt[dir]*dX[dir] + 0.5*dX[dir]; \
            } else {\
            }\
        } \
        proto_map_##name (a_x(0),a_x(1),a_x(2),X0[0],X0[1],X0[2],a_block); \
    } \
    const char* myname = #name; \
    inline Array<double, DIM> map(const Array<double, DIM>& X, unsigned int block) const { \
        Array<double, DIM> x; \
        proto_map_##name (x[0], x[1], x[2], X[0], X[1], X[2], block); \
        return x; \
    } \
    }; \
static struct_proto_map_##name name;

#elif DIM==2
#define PROTO_MAP_END(name) \
    } \
    struct struct_proto_map_##name { \
    template<MemType MEM = MEMTYPE_DEFAULT> \
    inline ACCEL_DECORATION \
    void operator() (Point& a_pt, Var<double, DIM, MEM>& a_x, \
            unsigned int a_block, \
            Point a_blockSize, \
            Centering a_centering) const \
    { \
        Array<double, DIM> dX = a_blockSize; \
        for (int dir = 0; dir < DIM; dir++) { dX[dir] = 1.0/dX[dir]; } \
        Array<double, DIM> X0; \
        for (int dir = 0; dir < DIM; dir++) \
        { \
            X0[dir] = a_pt[dir]*dX[dir]; \
            if (!((int)a_centering == dir || a_centering == PR_NODE)) \
            { \
                X0[dir] = X0[dir] + 0.5*dX[dir]; \
            } \
        } \
        double Z = 0; \
        proto_map_##name (a_x(0),a_x(1),Z,X0[0],X0[1],Z,a_block); \
    } \
    const char* myname = #name; \
    inline Array<double, DIM> map(const Array<double, DIM>& X, unsigned int block) const { \
        Array<double, DIM> x; \
        double Z = 0; \
        proto_map_##name (x[0], x[1], Z, X[0], X[1], Z, block); \
        return x; \
    } \
    }; \
static struct_proto_map_##name name;
#endif
template<typename Func, MemType MEM = MEMTYPE_DEFAULT>
class MBMap
{
    public:

    /// Constructor
    /**
      Create a MBMap and allocate space for data members

      \param a_layout  The MMB domain
      \param a_order   Order of accuracy for differentiation (default: 4)
      */
    inline MBMap(
            const Func& a_func,
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            Point                a_boundGhost,
            unsigned int a_order = 4);

    inline void define(
            const Func& a_func,
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            Point                a_boundGhost,
            unsigned int a_order = 4);

    /// Map Access
    inline MBLevelBoxData<double, DIM, MEM, PR_NODE>& map();
    
    /// Jacobian Access
    inline MBLevelBoxData<double, 1, MEM, PR_CELL>& jacobian();

    /// Compute Mapped Point From Array
    inline Array<double, DIM> operator()(
            const Array<double, DIM> a_X,
            unsigned int a_block) const;
    
    /// Compute Mapped Point From MBDataPoint
    inline Array<double, DIM> operator()(
            const MBDataPoint& a_pt,
            Point a_bound = Point::Zeros()) const;
    
    /// Compute Mapped Point From Point
    inline Array<double, DIM> operator()(
            const Point& a_cell,
            const Point& a_bound,
            unsigned int a_block) const;

    /// Compute Mapped Region
    inline BoxData<double, DIM> operator()(
            const Box& a_domain,
            unsigned int a_block,
            unsigned int a_refBlock,
            Centering a_ctr = PR_NODE) const;
    
    /// Compute Mapped Region (No Reference Block)
    inline BoxData<double, DIM> operator()(
            const Box& a_domain,
            unsigned int a_block,
            Centering a_ctr = PR_NODE) const;
    
    // Find Cell
    inline bool contains(
            const MBIndex& a_patch, const Point& a_cell, const Array<double, DIM>& a_X) const;
    private:

        unsigned int m_order;
        const Func* m_func;
        MBLevelBoxData<double, DIM, MEM, PR_NODE>   m_X;
        MBLevelBoxData<double, 1, MEM, PR_CELL>   m_J;
};

// TODO: This doesn't quite work without a copy function / move constructor
template<typename Func, MemType MEM = MEMTYPE_DEFAULT>
MBMap<Func, MEM>&& buildMap(
            const Func& a_func,
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            unsigned int a_order = 4)
{
    return MBMap<Func, MEM>(a_func, a_layout, a_ghost, a_order);
}

#else
    // Version 1

    template<typename T, MemType MEM>
    PROTO_KERNEL_START
    void f_determinantF(Var<T, 1, MEM>& a_det, Var<T, 3, MEM, 3>& a_mat)
    {
        a_det(0) = 0;
        for (int ii = 0; ii < 3; ii++)
        {
            a_det(0) += a_mat(ii%3, 0)*a_mat((1+ii)%3, 1)*a_mat((2+ii)%3, 2);
            a_det(0) -= a_mat(ii%3, 2)*a_mat((1+ii)%3, 1)*a_mat((2+ii)%3, 0);
        }
        //(+ 0,0  1,1  2,2) + (1,0  2,1  0,2) + (2,0  0,1  1,2)
        //(- 0,2  1,1  2,0) - (1,2  2,1  0,0) - (2,2  0,1  1,0)
    }
    PROTO_KERNEL_END(f_determinantF, f_determinant);

    template<typename T, MemType MEM = MEMTYPE_DEFAULT>
    class MBMap
    {
        public:
      
        /// Constructor
        /**
            Create a MBMap and allocate space for data members

            \param a_layout  The MMB domain
            \param a_order   Order of accuracy for differentiation (default: 4)
        */
        inline MBMap(MBDisjointBoxLayout& a_layout,
                std::array<Point, DIM+1>& a_ghost,
                unsigned int a_order = 4);
        inline void define(MBDisjointBoxLayout& a_layout,
                std::array<Point, DIM+1>& a_ghost,
                unsigned int a_order = 4);
        
        /// Compute Metric Terms
        /**
            Compute the map, cofactor matrices, and Jacobian from a user defined function
            a_func must have the signature:
            
            template<typename T, MemType MEM>
            void a_func(Point& a_pt, Var<T, 3, MEM>& a_X, unsigned int a_block, ...)

            Where <code>a_block</code> is the block index and ... are any optional arguments
            needed and which conform to the input requirements of Proto::forall. 
        */
        template<typename Func, typename... Srcs> 
        inline void compute(const Func& a_func, Srcs&&... a_srcs);

        /// Map Access
        inline MBLevelBoxData<T, 3, MEM, PR_NODE>& map();
        
        /// Jacobian Access
        inline const MBLevelBoxData<T, 1, MEM, PR_CELL>& jacobian() const;
        
        /// Compute Deformation Gradient
        inline BoxData<T, 3, MEM, 3> deformation(unsigned int a_blockIndex,
                const LevelIndex& a_patchIndex, int a_dir) const;
        inline BoxData<T, 3, MEM, 3> deformation(const MBIndex& a_index, int a_dir) const;
        
        /// Compute Cofactor Matrix
        inline BoxData<T, 3, MEM, 3> cofactors(unsigned int a_blockIndex,
                const LevelIndex& a_patchIndex, int a_dir) const;
        inline BoxData<T, 3, MEM, 3> cofactors(const MBIndex& a_index, int a_dir) const;
        
        // Find Cell
        inline bool contains(
                const MBIndex& a_patch, const Point& a_cell, const std::array<T, 3>& a_X) const;
        //inline std::vector<MBIndex> find(const std::array<T, 3>& a_X) const;
        //inline MBIndex find(const std::array<T, 3>& a_X, unsigned int a_blockIndex) const;
        private:
        
        unsigned int m_order;

        MBLevelBoxData<T, 3, MEM, PR_NODE>   m_X;
        MBLevelBoxData<T, 1, MEM, PR_CELL>   m_J;

        bool m_defined;  ///< Has data been allocated
        bool m_computed; ///< Has data been computed
    };
#endif // end Version 1
#include "implem/Proto_MBMapImplem.H"
} // end namespace Proto
#endif // end include guard
