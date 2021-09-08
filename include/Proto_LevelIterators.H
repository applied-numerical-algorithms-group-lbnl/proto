#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef Proto_LEVELITERATORS_H_
#define Proto_LEVELITERATORS_H_
#include <vector>
#include "Proto.H"
#include "Proto_DataIndex.H"
#include "Proto_DisjointBoxLayout.H"

namespace Proto
{
    /// Distributed Data Iterator
    /**
        Iterate through the local contents of a DisjointBoxLayout or LevelBoxData.
        Unline normal iterators, DataIterator dereferences into a DataIndex which can
        be used to index into a DisjointBoxLayout, LevelBoxData, or any other "level scoped"
        data holder.

        The intended syntax for for iterating is:

        DisjointBoxLayout layout(...);
        LevelBoxData<double> data(layout, ...);

        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            Box b = layout[*iter];
            BoxData<double>& patch = data[*iter];
            // etc
        }
    */
    class DataIterator
    {
        protected:
        
        DisjointBoxLayout m_dbl;
        DataIndex m_current;
        DataIndex m_end;
        
        public:
       
        /// Manual Constructor
        /**
            Construct a DataIterator from an existing DisjointBoxLayout.
            Generally there is no reason to do this; it's usually easier to call
            DisjointBoxLayout::begin() or LevelBoxData::begin(). 
            FIXME: Currently the iterator that is constructed is initialized as end() and I'm not sure why... -CLG
        */ 
        inline DataIterator(const DisjointBoxLayout& a_dbl);
        
        /// Set To Start
        /**
            Modifies the *this and returns it
        */
        inline DataIterator& begin();

        /// Set To End 
        /**
            modifies *this and returns it
        */
        inline DataIndex& end(){return m_end;};
        
        /// Continue Query
        /**
            Shortcut to check if *this has reached its end.
        */
        inline bool ok() const;

        /// Increment 
        inline DataIterator& operator++();
        
        /// Dereference
        /**
            Unlike normal iterators, dereferencing a DataIterator returns a DataIndex.
            DataIndex is accepted by all "level scope" data holders (e.g. LevelBoxData or DisjointBoxLayout).
        */
        inline const DataIndex& operator*() const;

        /// Get Patch Point
        /**
            Returns a Point associated with this patch in the level.
        */
        inline Point point() const;

        /// Returns the Box corresponding to that index. Undefined if invalid.
        /**
            Returns the Box associated with this patch in the layout.
            Does not take any possible LevelBoxData ghost cells into account. 
        */
        inline Box box() const; 
    };


    class NeighborIterator
    {
        protected:
        
        DisjointBoxLayout m_dbl;
        Box m_bx;
        BoxIterator m_bxit;
        Point m_shift;
        DataIndex m_current;
        DataIndex m_end;
        ProblemDomain m_dProbC;
        
        public:
        
        /// Constructor takes DBL, and Box that we are intersecting with the DBL.
        inline NeighborIterator(const DisjointBoxLayout& a_dbl, const Box& a_bx);
        //~NeighborIterator();

        /// Iterator functions begin, end, ++.
        inline NeighborIterator& begin();
        inline const DataIndex& end() const{return m_end;};
        inline bool ok() const {return (m_current != m_end);}
        inline NeighborIterator& operator++();

        /// Dereference operator returns the DataIndex corresponding to the current iterate.
        inline const DataIndex& operator*() const{return m_current;};
        /// The Point corresponding to the Box for that index. Undefined if invalid.
        inline Point point() const;
        /// The region to be copied from corresponding to that index. Undefined if invalid.
        inline Box srcBox() const;
        /// The region to be copied from corresponding to that index. Undefined if invalid.
        inline Box destBox() const;
        /// The amount by which the srcBox is shifted to obtain destBox a_bx.
        inline Point shift();
        // Processor ID for current iterate.
        inline int procID(){return m_dbl.procID(m_current);}; 
    };
#if 0
    /// Class that manages iteration over the flux register regions for a fine Box.
    /** 
      Used to to construct Copier for flux registers (Copier::FluxRegisterDefine).
      There is potentially one Coarse Box for each (LoHiSide,direction), that exists 
      only if the region defined by the register is not covered by another Box in the 
      fine DBL.
     */
    class FineFRIterator
    {
        public:
            /// Constructor takes DataIndex associated with the fine Box, and the refinement ratio.
            FineFRIterator(const DataIndex& a_fineIndex,
                    const Point& a_refRatio);
            ~FineFRIterator();

            /// Iterator functions begin, end, ++.
            inline void begin();
            inline DataIndex end();
            inline void operator++();
            /// Dereference operator returns the DataIndex corresponding to this object.
            inline DataIndex operator*{return m_fineIndex;};

            /// The low corner of the coarse Box corresponding to the current register.
            inline Point operator();

            /// The coarse Box corresponding to the current register.
            inline Box box();

            /// The shift applied to the coarse Box to get it into the valid region. 
            inline Point srcShift();
        protected:
            DataIndex m_fineIndex;
            pair<int,int> m_currentRegister;
            //first element is 0 or 1 (LoHiSide), second element is dir.
    };
    /// Class that manages iteration over the flux register regions for a coarse Box.
    class CoarseFRIterator
    {
        public:
            /// Coanstructor. 
            CoarseFRIterator(const& a_DisjointBoxLayout m_fineDBL,
                    const DataIndex& a_coarseDit,
                    const Point& a_refRatio);
            ~CoarseFRIterator();
            void begin();
            bool end();
            void operator++();
            DataIndex operator();
            Box box();
            Point srcShift();
            DataIndex fromIndex();
        protected:
            DataIndex m_coarseDit;
            DisjointBoxlayout m_fine;
            pair<LoHiSide,int> m_currentRegisterSet;
            BoxIterator m_bxit;
            // For each LoHiSide, direction, we iterate over a box corresponding to the
            // fine boxes coarsened to single Points inside the coarse box.
    }
#endif
#include "implem/Proto_LevelIteratorsImplem.H"
} // end Proto namespace.
#endif

