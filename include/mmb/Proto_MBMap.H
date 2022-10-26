#pragma once
#ifndef _PROTO_MB_MAP_
#define _PROTO_MB_MAP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"

namespace Proto
{

#if 1

#define PROTO_MAP_START(name) \
    ACCEL_DECORATION \
    void proto_map_##name (double& x, double& y, double& z, \
            const double& X, const double& Y, const double& Z, unsigned int block) {

#define PROTO_MAP_END(name) \
    } \
    struct struct_proto_map_##name { \
    template<MemType MEM = MEMTYPE_DEFAULT> \
    inline ACCEL_DECORATION \
    void operator() (Point& a_pt, Var<double, 3, MEM>& a_x, \
            unsigned int a_block, \
            Point a_blockSize) const \
    { \
        double X = (a_pt[0]*1.0)/(a_blockSize[0]*1.0); \
        double Y = (a_pt[1]*1.0)/(a_blockSize[1]*1.0); \
        double Z = (a_pt[2]*1.0)/(a_blockSize[2]*1.0); \
        proto_map_##name (a_x(0),a_x(1),a_x(2),X,Y,Z,a_block); \
    } \
    const char* myname = #name; \
    inline Array<double, 3> map(const Array<double, 3>& X, unsigned int block) const { \
        Array<double, 3> x; \
        proto_map_##name (x[0], x[1], x[2], X[0], X[1], X[2], block); \
        return x; \
    } \
    }; \
static struct_proto_map_##name name;

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
            unsigned int a_order = 4);

    inline void define(
            const Func& a_func,
            MBDisjointBoxLayout& a_layout,
            Array<Point, DIM+1>& a_ghost,
            unsigned int a_order = 4);

    /// Map Access
    inline MBLevelBoxData<double, 3, MEM, PR_NODE>& map();

    /// Compute Mapped Point
    inline Array<double, 3> operator()(const Array<double, 3> a_X, unsigned int a_block) const;

    // Find Cell
    inline bool contains(
            const MBIndex& a_patch, const Point& a_cell, const Array<double, 3>& a_X) const;
    private:

        unsigned int m_order;
        const Func* m_func;
        MBLevelBoxData<double, 3, MEM, PR_NODE>   m_X;
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
