#pragma once
#ifndef _PROTO_MB_MAP_
#define _PROTO_MB_MAP_

#include "Proto_MBDisjointBoxLayout.H"
#include "Proto_MBLevelBoxData.H"

namespace Proto
{

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

    template<typename T, MemType MEM>
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
        inline const MBLevelBoxData<T, 3, MEM, PR_NODE>& map() const;
        
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
#include "implem/Proto_MBMapImplem.H"
} // end namespace Proto
#endif // end include guard
