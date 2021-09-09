#pragma once
#ifndef _PROTO_AMR_DATA_
#define _PROTO_AMR_DATA_

#include "Proto_AMRGrid.H"
#include "Proto_LevelBoxData.H"

namespace Proto
{
    /// AMR Data Hierarchy
    /**
        A nested hierarchy of data defined on an AMRGrid.
    */
    template<typename T, unsigned int C=1, MemType MEM = MEMTYPE_DEFAULT, unsigned int CTR=0>
    class AMRData
    {
        public:
        
        /// Trivial Constructor
        AMRData() {m_defined = false;}

        /// Grid Constructor
        AMRData(AMRGrid& a_grid, Point a_ghost) { define (a_grid, a_ghost); }

        /// Single Grid Constructor
        /**
            Essentially a convenient overload of the AMRGrid constructor
        */
        AMRData(DisjointBoxLayout& a_layout, int a_maxLevels, Point a_ghost)
        {
            AMRGrid grid(a_layout, a_maxLevels);
            define(grid, a_ghost);
        }

        /// Define
        /**
            If *this is not yet defined, LevelBoxData are created based on the
            inputs. If *this was previously defined, redefine it based on the new
            AMRGrid (e.g. after AMRData::regrid is called).
        */
        inline void
        define(AMRGrid& a_grid, Point a_ghost);        
        
        /// Grid Access (Const)
        inline const AMRGrid& grid() const {return m_grid; }
        
        /// Grid Access (Non-Const)
        inline AMRGrid& grid() { return m_grid;}
        
        /// Level Data Access (Const)
        inline const LevelBoxData<T, C, MEM, CTR>&
        operator[](unsigned int a_level) const;
        
        /// Level Data Access (Non-Const)
        inline LevelBoxData<T, C, MEM, CTR>&
        operator[](unsigned int a_level);

        /// Zero Initialize
        inline void setToZero();

        /// Function Initialization
        /**
            Initializes *this using a function with a signature similar to that used by
            Proto::forallInPlace_p. This function is purely for convenience. 

            Func signature:
            void a_func(Point& a_point, Var<T, C, MEM>& a_thisData, double a_dx, ... )

            Unlike Proto::forall and its variants, this function cannot manage additional
            data holders (e.g. Srcs should not include a LevelBoxData or similar complex object).
            Unlike LevelBoxData::initialize(...), valid input functions must also include a 
            double argument representing the grid spacing, even if it is not used by the function
            itself. Elements of a_srcs should be plain-old-data such as single values,
            Points, or Boxes. 
        */
        template<typename Func, typename... Srcs>
        inline void initialize(double a_dx0, Func& a_func, Srcs... a_srcs);

        /// Recompute Grids
        /**
            Compute regridding tags using the data in *this and recompute
            the underlying AMRGrid using those tags. Then, redefine this 
            using the new AMRGrid.
        */
        inline void regrid(unsigned int a_level = 0);

        private:

        bool m_defined;
        AMRGrid m_grid;
        Point m_ghost;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        
    };
    
    typedef AMRData<char, 1, MEMTYPE_DEFAULT, 0> AMRTagData;

#include "implem/Proto_AMRDataImplem.H"
}
#endif //end include guard
