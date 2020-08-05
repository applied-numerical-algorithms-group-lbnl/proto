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
  class ProblemDomain
  {
  protected:
    Box m_bx;
    array<bool,DIM> m_isPeriodic;
    bool m_isDefined;
  public:
    ProblemDomain(){};
    //~ProblemDomain(){};
    ProblemDomain(const Box& a_bx,const array<bool,DIM> a_isPeriodic)
    {
      define(a_bx,a_isPeriodic);
    };
    inline void define(const Box& a_bx,const array<bool,DIM> a_isPeriodic)
    {
      m_isDefined = true;
      m_isPeriodic = a_isPeriodic;
      m_bx = a_bx;
      PROTO_ASSERT(m_bx.low() == Point::Zeros(),"Low corner of the problem domain not zero");
      // For the time being, we are assuming that the low corner of the problem domain is zero.
      // We can come back and fix this if necessary.
    };
                    
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
    inline Point shifted(const Point& a_pt) const
    {
      PR_assert(!((*this)&Box(a_pt,a_pt)).empty());
      return a_pt%(Point::Ones()*(m_bx.high() + 1));//Warning :: assumes low corner is set to zero
    };
    inline bool coarsenable(Point a_boxSize)
    {
      return m_bx.coarsenable(a_boxSize);
    };
    inline bool operator==(const ProblemDomain& a_input) const
    {
      return (m_bx== a_input.m_bx)&&(m_isPeriodic==a_input.m_isPeriodic)
              &&(m_isDefined==a_input.m_isDefined);
    };
    inline Box box() const {return m_bx;};
    inline ProblemDomain coarsen(Point a_refRatio) const
    {
      ProblemDomain retval;
      retval.m_bx = m_bx.coarsen(a_refRatio);
      retval.m_isPeriodic = m_isPeriodic;
      retval.m_isDefined = m_isDefined;
      return retval;
    };
    inline ProblemDomain refine(Point a_refRatio) const
    {
      ProblemDomain retval;
      retval.m_bx = m_bx.refine(a_refRatio);
      retval.m_isPeriodic = m_isPeriodic;
      retval.m_isDefined = m_isDefined;
      return retval;
    };
  };
  // end Proto namespace.
}  
#endif
