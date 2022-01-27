#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifndef _PROTO_PROBLEMDOMAIN_H_
#define _PROTO_PROBLEMDOMAIN_H_
#include <array>
#include "Proto.H"

using namespace std;
namespace Proto
{
    /// Represents a rectangular domain over which a problem can be defined, including periodic images.
    /**
        A ProblemDomain is specified by a Box defining the rectangular set of Points over which a 
        problem may be defined, plus an array<bool,DIM> specifying whether or not each of the
        coordinate directions is periodic. The Box is assumed to have the low corner to be all Zeros,
        and specifies the overall range of indices that contain all the boxes in a DisjointBoxLayout. 
    */
    class ProblemDomain
    {
        protected:
        
        Box m_bx;
        array<bool,DIM> m_isPeriodic;
        bool m_isDefined;
        
        public:
        
        ProblemDomain(){};
        
        /// Anisotropic Constructor
        /**
            Constructs a ProblemDomain with Box a_bx, and periodic directions given by a_periodic.
        */
        inline ProblemDomain(const Box& a_bx,const array<bool,DIM> a_isPeriodic);
        
        /// Isotropic Constructor
        inline ProblemDomain(const Box& a_bx, bool a_isPeriodic);
        
        /// Define
        inline void define(const Box& a_bx,const array<bool,DIM> a_isPeriodic);

        /// Box Intersection
        /**
            Returns the Box that is the intersection of a_bx with m_bx
            and all of its periodic images. Points in a_bx which have no periodic
            image are truncated.
        */
        inline Box operator&(Box a_bx) const;
       
        /// Periodic Image
        /**
            Converts input point into a periodic image IF such an image exists.
            Returns true if a_pt is inside this or a periodic image of this. Returns false
            otherwise.
        */ 
        inline bool image(Point& a_pt) const;
       
        /// Periodic Image (Original)
        /**
            This function does some of what ProblemDomain::image does. The main difference is
            that this function will throw an error if a_pt is outside of this domain
            or any of it's periodic images.
        */ 
        inline Point shifted(const Point& a_pt) const;
        
        /// Check If Coarsenable
        /**
            Must return true for a DisjointBoxLayout defined over *this to use fixed-sized
            boxes of size a_boxsize.
        */
        inline bool coarsenable(Point a_boxSize) const
        {
            return m_bx.coarsenable(a_boxSize);
        };

        /// Query Size
        /**
            Size in number of points in each direction.
        */
        inline Point size() const
        {
            return m_bx.high()  - m_bx.low() + Point::Ones();
        }
        
        /// Equality
        /**
            Two ProblemDomain objects are "==" if they have the same
            Box domain and periodicity.
        */
        inline bool operator==(const ProblemDomain& a_input) const
        {
            return (m_bx== a_input.m_bx)&&(m_isPeriodic==a_input.m_isPeriodic)
                &&(m_isDefined==a_input.m_isDefined);
        };

        /// Get Box
        inline Box box() const {return m_bx;};
        
        /// Get Periodicity
        inline array<bool,DIM> periodicflags() const{return m_isPeriodic;};
       
        inline bool isPeriodic(int a_dir) {return m_isPeriodic[a_dir]; }

        /// Coarsen
        /**
            Returns a ProblemDomain with a box given by this->box.coarsen(a_refRatio).
            Fails if this->coarsenable(a_refRatio) is false
         */
        inline ProblemDomain coarsen(Point a_refRatio) const;
        inline ProblemDomain coarsen(int a_refRatio) const;

        /// Refine
        /**
            Returns a ProblemDomain with a box given by this->box.refine(a_refRatio).
        */
        inline ProblemDomain refine(Point a_refRatio) const;
        inline ProblemDomain refine(int a_refRatio) const;
    }; // end class ProblemDomain

    /// Stream output for ProblemDomain.
    inline std::ostream& operator<< (std::ostream& os, const ProblemDomain& a_pd)
    {
        os << "ProblemDomain: Box = " << a_pd.box() <<
            " , periodicity = [" ;
        for (int d = 0; d < DIM; d++)
        {
            os<< a_pd.periodicflags()[d] ;
            if (d < DIM-1) os << " , ";
        }
        os << "]";
        return os;
    }
#include "implem/Proto_ProblemDomainImplem.H"
} // end namespace Proto  
#endif //end include guard
