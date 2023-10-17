
#ifndef __Proto_Transition_Utilities_H__
#define __Proto_Transition_Utilities_H__
#include <vector>
#include "../base/Proto_MayDay.H"
#include "../base/Proto_DisjointBoxLayout.H"
///
/**
   This is a set of utilities to allow users of other infrastructures
   to join the fun.
   This class has no data and is just standalone static functions.
   I am probably using way too many shared pointers;
   (some things have default constructors, others do not).
*/
namespace Proto_Transition
{
  ///
  /**
    This class has no data--just static functions.
   **/
  class MetaUtilities
  {
  public:

    typedef  Proto::MayDay<void>                     pr_mayday;
    typedef  Proto::Point                            pr_pt;
    typedef  Proto::Box                              pr_box;
    typedef  Proto::DisjointBoxLayout                pr_dbl;
    typedef  Proto::ProblemDomain                    pr_dom;

    ///Get a Proto::DisjointBoxLayout that conforms to the input meta data
    /** 
        A box here is its point refined by the box size.   
        pair.first = the valid box of the layout.
        pair.second = the processor assignment of said box.
        The input layout can be in any order or have any distribution but 
        the boxes MUST BE DISJOINT.
    **/
    static shared_ptr<pr_dbl> getLayout(
      const std::vector<std::pair<pr_pt, int> > &  a_layout,
      const pr_pt                               &  a_gridSize,
      const pr_dom                              &  a_domain)
    {
      /**
         Ignoring processor arguments and the fact that order can change.
         This will very much fail in parallel and will probably fail in serial.
      **/
      vector<pr_pt> coarsenedPatches(a_layout.size());
      for(int ibox  = 0; ibox < a_layout.size(); ibox ++)
      {
        coarsenedPatches[ibox] = a_layout[ibox].first;
      }
      shared_ptr<pr_dbl> retval(new pr_dbl(a_domain, coarsenedPatches, a_gridSize));
      return retval;
    }
  };

  ///
  /**
    This class has no data--just static functions.
   **/
  template <int ncomp>
  class DataUtilities
  {
  public:

    typedef  Proto::MayDay<void>                     pr_mayday;
    typedef  Proto::Point                            pr_pt;
    typedef  Proto::Box                              pr_box;
    typedef  Proto::DisjointBoxLayout                pr_dbl;
    typedef  Proto::ProblemDomain                    pr_dom;
    typedef  Proto::LevelBoxData<double, ncomp>      pr_level_data;
    typedef  Proto::BoxData<     double, ncomp>      pr_box_data;


    ///
    /**
     *  \param a_host_data - one buffer per box
     **/ 
    inline void
    copyToHost(
      std::vector< double * >  &    a_host_data,
      const pr_level_data      &    a_devi_data)
    {
      pr_mayday::Error("Doh! copyToHost not implemented");
    }

    ///
    /**
     *  \param a_devi_data  device data to be filled in
     *  \param a_host_data  matching host pointers for each BoxData in a_devi_data
     *  The host data has to be in fortran ordering and must contain the same footprint 
     *  as the device BoxData (it has to include ghost cell memory).
     **/
    inline void
    copyToDevice(
      pr_level_data                  &    a_devi_data,
      const std::vector< double * >  &    a_host_data)
    {
      pr_mayday::Error("Doh! copyToDevice not implemented");
    }

    ///
    /**
     *      \a_domain - problem domain
     *      \a_domain_face_dir is the direction of the domain face
     *      \a_lo_or_hi_side  = 0 for low side of the domain, = 1 for high.
     **/
    inline void
    domainGhostReflectEven(
      pr_level_data             &    a_devi_data,
      const pr_box              &    a_domain,
      int a_domain_face_dir,  int    a_lo_or_hi_side)  
    {
      pr_mayday::Error("Doh!  ghostReflectEven not implemented");
    }

    ///
    /**
     * \param     \a_domain - problem domain
     * \param     \a_domain_face_dir is the direction of the domain face
     * \param     \a_lo_or_hi_side  = 0 for low side of the domain, = 1 for high.
     **/
    inline void
    domainGhostReflectOdd(
      pr_level_data             &    a_devi_data,
      const pr_box              &    a_domain,
      int a_domain_face_dir,  int lo_or_hi_side)  
    {
      pr_mayday::Error("Doh!  ghostReflectOdd not implemented");
    }

    ///
    /**
     *  \param a_fine_data - data whose ghost is getting filled
     *  \param a_coar_data - data   at the next coarser level
     *  \param a_domain    - domain at the next coarser level
     *  \param ref_ratio   - refrat to the next coarser level
     **/
    template <int interpolation_order>
    inline void
    fillGhostOverCoarseFineBoundaries(
      pr_level_data             &    a_fine_data,
      const pr_level_data       &    a_coar_data,
      const pr_box              &    a_domain,
      const pr_pt               &    a_ref_ratio)
    {
      pr_mayday::Error("Doh!  fillGhostOverCoarseFineBoundaries not implemented");
    }
    
  };
}
#endif
