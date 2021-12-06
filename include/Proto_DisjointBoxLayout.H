#pragma once

#ifndef _PROTO_DisjointBoxLayout_H_
#define _PROTO_DisjointBoxLayout_H_

#include "Proto_MayDay.H"
#include "Proto_Point.H"
#include "Proto_Box.H"
#include "implem/Proto_Stack.H"
#include "Proto_ProblemDomain.H"
#include "Proto_DBLInternals.H"
#include "Proto_DataIndex.H"
#include "Proto_Morton.H"
#include "Proto_SPMD.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <iomanip>
#include <stack>
#include <memory>
#include <array>
#include <unordered_map>

//using namespace std;
//using std::shared_ptr;
//using std::array;
//using std::vector;
//using std::stack;
//using namespace Proto;

namespace Proto 
{ 
    // Forward Declarations
    class DataIterator; 

    /// Disjoint Box Layout
    /**
        Layout of disjoint, logically rectangular regions of DIM-dimensional space.
        This object is a simplified generalization of the DisjointBoxLayout class
        used in Chombo, now specialized to use fixed box sizes. 
    */
  class DisjointBoxLayout
    { 
      friend class DataIndex;
        
    private:
      ProblemDomain             m_problemDomain;
      Point                     m_boxSize;  
      DataIndex                 m_end;     
      shared_ptr<DBLInternals>  m_internals;
        /* 
           The use of  shared_ptr<DBLInternals> makes the proliferation of copies not be a 
           performance barrier, as well as making the data pointed at by m_internals reusable 
           under coarsening or refinement.If making m_problemDomain and m_boxSize objects 
           becomes a performance issue, can wrap them in a shared pointer as well.
         */
        
        public:
        
        /// 
        inline DisjointBoxLayout();

        /// Constructor (Full Domain)
        /**
            Creates a DisjointBoxLayout by covering an input ProblemDomain completely with boxes of
            a predetermined size. If MPI is used, boxes will be distributed over available
            processors using Morton ordering.
        */
        inline DisjointBoxLayout(const ProblemDomain   & a_problemDomain, 
                const Point           & a_boxSize);

        /// Constructor (Partial Domain)
        /**
            Creates a DisjointBoxLayout by partially covering an input ProblemDomain with boxes
            according to a set of Points.

            Each element of a_coarsenedPatches represents a box in the layout which has been coarsened
            by boxSize (resulting in a Box with a single element). 
           
            TODO: Needs LaTeX writeup 
            Concretely, for each Point p in the input array, a Box b will be constructed in the layout
            such that:
            
            Point p => Box b = Box(p,p).refine(boxSize)
                             = [boxSize*p, boxSize*(p + (1,1,...)) - (1,1,...)]
            
            \param a_problemDomain   ProblemDomain containing the DBL.
            \param a_coarsenedPatches   contains the Points correcpoinding the patches in the DBL for all the processors.
            \param a_boxSize is the size of a patch. Note that this allows anisotropic patches.
        */
        inline DisjointBoxLayout(const ProblemDomain   & a_problemDomain,
                const vector<Point>   & a_coarsenedPatches,
                const Point           & a_boxSize);

        /// Copy Constructor
        /**
            Copies share pointers to the same internal structure
        */
        inline DisjointBoxLayout(const DisjointBoxLayout& a_input)
        {
            m_internals = a_input.m_internals;
            m_problemDomain = a_input.m_problemDomain;
            m_boxSize = a_input.m_boxSize;
            m_end = a_input.m_end;
        };
        
        /// Define (Full Domain)
        /**
            Construct a DisjointBoxLayout lazily or redefine an existing one. 
        */
        inline void define(const ProblemDomain   & a_problemDomain, 
                           const Point           & a_boxSize);
        
        ///  Define (Partial Domain)
        /**
            Construct a DisjointBoxLayout lazily or redefine an existing one. 
        */
        inline void define(const ProblemDomain   & a_problemDomain,
                const vector<Point>   & a_coarsenedPatches,
                const Point           & a_boxSize);
        
        /// Assignment
        /**
            Copies share pointers to the same internal structure
        */
        inline DisjointBoxLayout& operator=(const DisjointBoxLayout& a_input)
        {
            m_internals = a_input.m_internals;
            m_problemDomain = a_input.m_problemDomain;
            m_boxSize = a_input.m_boxSize;
            m_end = a_input.m_end;
            return *this;
        };
        
        /// Equality Operator
        inline bool operator==(const DisjointBoxLayout& a_input) const
        {
            return ((a_input.m_internals == m_internals)&&
                    (a_input.m_problemDomain == m_problemDomain)&&
                    (a_input.m_boxSize == m_boxSize));
        };

        /** @name Access operators. 
         */

        /// Get Internals
        /**
            TODO: try to solve this with friends -CLG
            Needed by DataIndex -PC
        */
        inline const shared_ptr<DBLInternals> dblInternals() const{return m_internals;};

        /// Index Access
        /**
            Return the Box associated with a DataIndex by calling this->box(a_index)
        */
        inline Box operator[](const DataIndex& a_index) const;

        /// Box Indexing
        /**
            Return the Box associated with a DataIndex.
        */
        inline Box box(const DataIndex& a_index) const;
        
        /// Point Indexing
        /**
            Return the Point in patch space associated with a DataIndex
        */
        inline Point point(const DataIndex& a_index) const;

        /// Get Process ID
        inline int procID(const DataIndex& a_index) const; //consistent syntax
        inline int procid(const DataIndex& a_index) const;

        /// Patch Offset
        /**
            Returns the linear offset (in number of boxes) of the first patch on proc a_proc
        */
        inline unsigned int offset(int a_proc) const;
        inline unsigned int offset() const { return offset(Proto::procID()); }

        /// Get Index From Point
        /**
            Find the DataIndex associated with an input patch point.
            This function is the inverse of DisjointBoxLayout::point(DataIndex).

            Points outside the ProblemDomain will still correspond to a DataIndex provided
            they are in periodic images of the ProblemDomain.
        */
        inline DataIndex index(Point& a_coarsenedPt) const;

        /// Size
        /**
            Return the number of Boxes contained in *this including those on other processors.
        */
        inline unsigned int size() const;

        /// Local Size
        /**
            Return the number of Boxes on this processor.
        */
        inline unsigned int localSize() const;
        inline unsigned int mySize() const;

        /// Global Indexing
        /**
            Get the DataIndex associated with a global integer index.
        */
        inline DataIndex index(unsigned int a_indexInt) const;

        /// Local Indexing
        /**
            Get the DataIndex associated with a local integer index.
        */
        inline DataIndex localIndex(unsigned int a_myIndexInt) const;
        inline DataIndex myIndex(unsigned int a_myIndexInt) const;

        /// Local Indexing (Inverse)
        /**
            Given a DataIndex, compute the associated local integer index.
            Results in an error if the DataIndex is out of bounds for
            this processor. 
        */
        inline int localIntIndex(const DataIndex& a_di) const;
        inline int myIntIndex(const DataIndex& a_di) const;

        /// Get Box Size
        /**
            Return the constant size of all boxes in the layout.
        */
        inline Point boxSize() const{return m_boxSize;};

        /// Get Problem Domain
        inline ProblemDomain domain() const {return m_problemDomain;};
        inline ProblemDomain problemDomain() const{return m_problemDomain;};

        /// Get Problem Domain in Patch Space
        /**
            Returns the problem domain box in patch-space.
        */
        inline ProblemDomain patchDomain() const {return domain().coarsen(boxSize()); }

        /// Get Max Local Size
        /**
            Returns the maximum size, in number of patches, on any processor.
            TODO: Currently unused; might be useful for load balancing.
        */
        inline int maxLocalSize() const;
        
        /// Get Min Local Size
        /**
            Returns the minimum size, in number of patches, on any processor.
            TODO: Currently unused; might be useful for load balancing.
        */
        inline int minLocalSize() const;

        /// Find Index
        /**
            Given a patch-point, return the DataIndex of the associated Box if
            it is in the layout.
            TODO: Someone check me on this description -CLG
        */
        
        inline DataIndex find(const Point& a_pt) const;
        ///@}
        /// Iterator End
        inline const DataIndex& end() const
        {
            return m_end;
        };
        
        /// Get Iterator
        inline DataIterator begin() const;
        
        /// Coarsenable Query
        /**
            Checks if the layout is coarsenable by a (possibly anisotropic) refinement ratio.
            This property requires the ProblemDomain to be coarsenable as well as
            boxSize() % refRatio == (0,0,....)
        */
        inline bool coarsenable(const Point& a_refrat) const;
        
        /// Simple Coarsen
        /**
            Coarsens the layout by an input (possibly anisotropic) refinement ratio.
            coarsenable(refRatio) must be true or this function will result in an error. 

            The coarsened layout has:
            coarseLayout.problemDomain().box() == fineLayout.problemDomain.box().coarsen(refRatio)
            coarseLayout[ii] == fineLayout[ii].coarsen(refRatio)   for each ii in [0, size())
            TODO: Write this in LaTeX
            TODO: Can this change the number of Boxes? -CLG
        */
        inline DisjointBoxLayout coarsen(const Point& a_refrat) const;
        
        /// Coarsen
        /**
            Creates a new DisjointBoxLayout with a specified box size that has been coarsened by 
            a valid ratio. coarsenable(refRatio) must be true or this function will result in an error.
        */
        inline DisjointBoxLayout coarsen(const Point& a_refrat, const Point& a_boxSize) const;

        /// Refine
        /**
            Refines the layout by an input (possibly anisotropic) refinement ratio.
            Unline coarsen(...), this function should always be possible.

            The refined layout has:
            fineLayout.problemDomain().box() == coarseLayout.problemDomain.box().refine(refRatio)
            fineLayout[ii] == coarseLayout[ii].refine(refRatio)   for each ii in [0, size())
            TODO: Write this in LaTeX
            TODO: If there is a max boxSize, this could change the number of Boxes -CLG
        */
        inline DisjointBoxLayout refine(const Point& a_refrat) const;

        /// Check Radial Symmetry
        /**
            For debugging mesh refinement procedures
        */
        inline bool radialSymmetry(Point a_origin) const;
        
        /// Check Mirror Symmetry
        /**
            For debugging mesh refinement procedures
        */
        inline bool mirrorSymmetry(Point a_origin, int a_coord) const;

        /// Print
        inline void print(std::string a_name = "") const;

        template <MemType MEMTYPE>
        static StackAlloc<MEMTYPE>& getStackAllocator()
        {
            static StackAlloc<MEMTYPE> s_allocator;
            return s_allocator;
        }
#ifdef PROTO_CUDA
        //Returns the single stream
        static protoStream_t getCurrentStream()
        {
            static protoStream_t currStream;
            static bool init=false;
            if(!init)
            {
                protoStreamCreate(&currStream);
                init=true;
            }
            return currStream;
        }
#endif
};

/// Stream output for DBL.
inline std::ostream& operator<< (std::ostream& os, const DisjointBoxLayout& a_dbl)
{
    os << "DisjointBoxLayout: \n" ;
    os << a_dbl.problemDomain() << "\n";
    auto internals = a_dbl.dblInternals();
    os << "coarsened Domain = " << internals->coarsenedDomain << "\n" ;
    os << "Points and ranks: \n";
    for (int k = 0; k < internals->allBoxes.size();k++)
    {
        os << internals->allBoxes[k].first 
            << " , " << internals->allBoxes[k].second << "\n";
    }
    os << "StartProc vector: \n";
    for (int k = 0; k < internals->startProc.size() ; k++)
    {
        os << "Starting patch number for processor " << k 
            << " = " << internals->startProc[k] << "\n";
    }
    return os;
}
#include "implem/Proto_DisjointBoxLayoutImplem.H"
} //end namespace Proto

#endif
