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
      A ProblemDomain is specified by a Box defining the rectangular set of Points over which a problem may be defined, plus an array<bool,DIM> specifying whether or not each of the coordinate directions is periodic. The Box is assumed to have the low corner to be all Zeros, and specifies the overall range of indices that contain all the boxes in a DisjointBoxLayout. 
     */
    class ProblemDomain
    {
        protected:
            Box m_bx;
            array<bool,DIM> m_isPeriodic;
            bool m_isDefined;
        public:
            ProblemDomain(){};
            //~ProblemDomain(){};
            
            /// Constructor constructs a ProblemDomain with Box a_bx, and periodic directions given by a_periodic.
            ProblemDomain(const Box& a_bx,const array<bool,DIM> a_isPeriodic)
            {
                //FIXME
                PROTO_ASSERT(a_bx.low() == Point::Zeros(),
                    "ProblemDomain constructor | Error: Many functions in ProblemDomain do not work correctly if the low corner is not (0,...,0).");
                define(a_bx,a_isPeriodic);
            };

            /// Isotropic Constructor
            ProblemDomain(const Box& a_bx, bool a_isPeriodic)
            {
                std::array<bool, DIM> periodic;
                for (int ii = 0; ii < DIM; ii++){periodic[ii] = a_isPeriodic;}
                define(a_bx,periodic);
            };
            
            /// Define - same arguments as constructor, but allows weak construction.
            inline void define(const Box& a_bx,const array<bool,DIM> a_isPeriodic)
            {
                m_isDefined = true;
                m_isPeriodic = a_isPeriodic;
                m_bx = a_bx;
                PROTO_ASSERT(m_bx.low() == Point::Zeros(),"Low corner of the problem domain not zero");
                // For the time being, we are assuming that the low corner of the problem domain is zero.
                // We can come back and fix this if necessary.
            };

            /// Intersection of problem domain with a Box a_bx, including periodicity. 
            /** Returns the Box that is the intersection of a_bx with m_bx and all of its periodic images. */
            inline Box operator&(Box a_bx) const
            {
                Point lowCorner = a_bx.low();
                Point highCorner = a_bx.high();
                for (int dir = 0; dir < DIM; dir++)
                {
                    if (!m_isPeriodic[dir])
                    {
                        lowCorner[dir] = max(lowCorner[dir],m_bx.low()[dir]);
                        highCorner[dir] =  min(highCorner[dir],m_bx.high()[dir]);
                    }
                }
                return Box(lowCorner,highCorner);
            };
           
            /// Periodic Image
            /**
                Converts input point into a periodic image IF such an image exists.
            */ 
            inline bool image(Point& a_pt) const
            {
                if (m_bx.empty()) { return false; }
                if (m_bx.contains(a_pt)) { return true; }
                
                Box domain = m_bx.shift(-m_bx.low());
                Point point = a_pt - m_bx.low();
                for (int ii = 0; ii < DIM; ii++)
                {
                    if (point[ii] < domain.low()[ii] || point[ii] > domain.high()[ii])
                    {
                        if (!m_isPeriodic[ii]) { return false; }
                        else {
                            int a = point[ii]; int b = domain.size(ii);
                            int pi = a % b;
                            if (pi < 0) { pi += m_bx.size(ii); }
                            point[ii] = pi;
                        }
                    }
                }
                a_pt = point + m_bx.low();
                return true;
            }
            
            inline Point shifted(const Point& a_pt) const
            {
                PR_assert(  !( (*this) & Box(a_pt,a_pt)  ).empty() );
                return a_pt % (Point::Ones()*(m_bx.high() + 1));
                //FIXME: Assumes low corner is (0,...,0). Also doesn't check which directions are periodic. 
            };
            /// Checks to see whether the ProblemDomain Box is coarsenable by a_boxsize.
            /**
              Must return true for a DisjointBoxLayout defined over *this to use fixed-sized boxes of size a_boxsize.
             */
            inline bool coarsenable(Point a_boxSize) const
            {
                return m_bx.coarsenable(a_boxSize);
            };
            inline Point size() const{return m_bx.high()+Point::Ones();}; //FIXME: This only works of low = (0,...,0)
            inline bool operator==(const ProblemDomain& a_input) const
            {
                return (m_bx== a_input.m_bx)&&(m_isPeriodic==a_input.m_isPeriodic)
                    &&(m_isDefined==a_input.m_isDefined);
            };
            /// Returns Box that defines the *this.
            inline Box box() const {return m_bx;};
            /// Returns the array of bools defining the periodic directions.
            inline array<bool,DIM> periodicflags() const{return m_isPeriodic;};
            /// Returns a ProblemDomain with a box given by this->box.coarsen(a_refRatio).
            /** this->coarsenable(a_refRatio) must be true, and is enforced by a PR_assert.
             */
            inline ProblemDomain coarsen(Point a_refRatio) const
            {
                PR_assert(this->coarsenable(a_refRatio));
                ProblemDomain retval;
                retval.m_bx = m_bx.coarsen(a_refRatio);
                retval.m_isPeriodic = m_isPeriodic;
                retval.m_isDefined = m_isDefined;
                PR_assert(m_isDefined);
                return retval;
            };
            /// Returns a ProblemDomain with a Box given by this->box.coarsen(a_refRatio).
            inline ProblemDomain refine(Point a_refRatio) const
            {
                ProblemDomain retval;
                retval.m_bx = m_bx.refine(a_refRatio);
                retval.m_isPeriodic = m_isPeriodic;
                retval.m_isDefined = m_isDefined;
                return retval;
            }; 
    };
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
    // end Proto namespace.
}  
#endif
