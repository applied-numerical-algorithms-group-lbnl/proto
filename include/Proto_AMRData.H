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
    template<typename T, unsigned int C=1, MemType MEM = MEMTYPE_DEFAULT, Centering CTR = CELL>
    class AMRData
    {
        public:
        
        /// Trivial Constructor
        AMRData() {m_defined = false;}

        /// Grid Constructor
        AMRData(AMRGrid& a_grid, Point a_ghost) { define (a_grid, a_ghost); }

        /// Single Grid Constructor
        /**
            A convenient overload. Creates an AMRGrid with a single level based on a_layout.
        */
        AMRData(DisjointBoxLayout& a_layout, int a_maxLevels, Point a_ghost)
        {
            AMRGrid grid(a_layout, a_maxLevels);
            define(grid, a_ghost);
        }

        /// Define
        /**
            If *this is not yet defined, LevelBoxData are created based on the
            inputs.
            TODO:
            If *this was previously defined, redefine it based on the new
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
        
        /// Function Initialization (With Convolution)
        /**
            Initialize *this using a user defined function for the cell-centered quantity
            and then convert the result to a cell-averaged quantity through 4th order
            convolution. See <code>AMRData::initialize</code> for further documentation.
        */
        template<typename Func, typename... Srcs>
        inline void initConvolve(double a_dx0, Func& a_func, Srcs... a_srcs);

        /// Average Down Fine Levels
        /**
            Synchronizes data between levels by recursively replacing coarse level data with the geometric
            average of overlying fine data where it exists.

            TODO: figure out how this should be defined for CTR != 0 (if at all)
        */
        inline void averageDown();
        
        /// Interpolate Fine Boundaries
        /**
            Populate coarse-fine boundary ghost cells with data interpolated from the next 
            coarser level for all levels > 0.
             
            TODO: figure out how this should be defined for CTR != 0 (if at all)
        
            \param a_order  Order of interpolation. Valid values: 3, 5 (default 5)
        */
        inline void interpBoundaries(int a_order = 5);

        /// Compute Integral
        /**
            Compute the discrete integral over all valid data.

            \param a_cdx    Coarsest grid spacing (Level 0)
        */
        inline double integrate(double a_cdx, unsigned int a_c = 0) const;
        
        /// Compute Max Norm
        /**
            Compute the max norm of all valid data
            
            \param a_cdx    Coarsest grid spacing (Level 0)
        */
        inline double absMax(unsigned int a_c = 0) const;
       
        /// Exchange
        /**
            Calls LevelBoxData::exchange for each level
        */
        inline void exchange();
            
        /// Recompute Grids
        /**
            TODO:
            Compute regridding tags using the data in *this and recompute
            the underlying AMRGrid using those tags. Then, redefine this 
            using the new AMRGrid.

            a_level specifies the base level; only data / grids finer than this level
            will be affected, and proper nesting will not be enforced on a_level or any
            coarser levels. 
        */
        inline void regrid(unsigned int a_level = 0);

        /// Compute Tags
        /**
            Populate tags using the magnitude of the undivided difference of the data.
            This is only public for testing and will later be called as a subroutine
            of AMRData::regrid(...)
            TODO:
            Specify check all components instead of just the 0 component.
            Use divided difference / mean value
        */
        inline void computeTags(LevelTagData& a_tags,
            unsigned int a_level,
            Point a_bufferSize,
            T a_threshold) const; 

        /// Get Number of Levels
        /**
            Returns the number of levels in *this with a defined LevelBoxData.
            This value will in general be different from the maximum number of levels.
            TODO:
            if numLevels() and grid().numLevels() are different, this is probably a bug. -CLG
        */
        inline int numLevels() const { return m_data.size(); }

        /// Get Max Number of Levels
        inline int maxLevels() const { return m_grid.maxLevels(); }

        /// 
        private:

        inline void buildInterpStencils(int a_order);

        bool m_defined;
        AMRGrid m_grid;
        Point m_ghost;
        std::vector<std::shared_ptr<LevelBoxData<T, C, MEM, CTR>>> m_data;
        
        InterpStencil<double> m_boundInterp; 
        int m_boundInterpOrder;
    };
    
    typedef AMRData<char, 1, MEMTYPE_DEFAULT, CELL> AMRTagData;

#include "implem/Proto_AMRDataImplem.H"
}
#endif //end include guard
