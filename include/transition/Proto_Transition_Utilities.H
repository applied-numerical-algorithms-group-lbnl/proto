
#ifndef __Proto_Transition_Utilities.H__
#define __Proto_Transition_Utilities.H__

#include "Proto.H"
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
  using std::shared_ptr;
  using std::shared_vector;
  using std::pair;
  typedef  Proto::MayDay                           pr_mayday;
  typedef  Proto::Point                            pr_point;
  typedef  Proto::Box                              pr_box;
  typedef  Proto::DisjointBoxLayout                pr_dbl;
  template <int ncomp>
  class DataUtilities
  {
  public:
    typedef  Proto::LevelBoxData<double, ncomp>      pr_level_data;
    typedef  Proto::BoxData<     double, ncomp>      pr_box_data;


    ///
    /** get me a dbl based on box - proc pairs
     **/
    inline pr_dbl getLayout(
            const vector<pair<pr_box, int> > >&  a_layout)
    {
      shared_ptr<pr_dbl> retval;
      PROTO_ASSERT(a_layout.size() > 0,
              "DataUtilities::getLayout | Error: empty input");
      
      auto sizes = a_layout[0].first.sizes();
      vector<pr_point> patches;
      for (auto iter : a_layout)
      {
        // the following assertion can be relaxed if necessary, but it simplifies the math without losing much
        PROTO_ASSERT(iter.first.low().min() >= 0,
                "DataUtilities::getLayout | Error: input boxes must be in the strictly non-negative index space");
        // check for constant size / shape
        PROTO_ASSERT(iter.first.sizes() == sizes,
                "DataUtilities::getLayout | Error: inputs do not represent a constant box size");
        // check for overlapping / improperly shifted boxes
        PROTO_ASSERT(iter.first.low() % sizes == Point::Zeros(),
                "DataUtilities::getLayout | Error: detected overlapping boxes");
        patches.push_back(iter.first.low() / sizes);
      }


      pr_mayday::Error("Doh! getMetaData not implemented")
      return retval;
    }

    ///
    /**
     *  \param a_host_data - one buffer per box
     **/ 
    inline void
    copyToHost(
            vector< double * >  &    a_host_data,
            const pr_level_data &    a_devi_data)
    {
      pr_mayday::Error("Doh! copyToHost not implemented")
    }

    ///
    /**
     *  \param a_devi_data - one buffer per box
     **/
    inline void
    copyToDevice(
            pr_level_data             &    a_devi_data,
            const vector< double * >  &    a_host_data)
    {
      pr_mayday::Error("Doh! copyToDevice not implemented")
    }

    ///
    /**
     *      \a_domain - problem domain
     **/
    inline void
    domainGhostReflectEven(
            pr_level_data             &    a_devi_data,
            const pr_box              &    a_domain)
    {
      pr_mayday::Error("Doh!  ghostReflectEven not implemented");
    }

    ///
    /**
     *      \a_domain - problem domain
     **/
    inline void
    domainGhostReflectOdd(
            pr_level_data             &    a_devi_data,
            const pr_box              &    a_domain)
    {
      pr_mayday::Error("Doh!  ghostReflectOdd not implemented");
    }

    ///
    /**
     *  \param pr_coar_box - coarse problem domain
     *  \param ref_ratio
     **/
    template <int interpolation_order>
    inline void
    fillGhostOverCoarseFineBoundaries(
            pr_level_data             &    a_fine_data,
            const pr_level_data       &    a_coar_data,
            const pr_box              &    a_domain,
            const pr_point            &    a_ref_ratio)
    {
      pr_mayday::Error("Doh!  fillGhostOverCoarseFineBoundaries not implemented");
    }
    
  };
}
#endif
