
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
  typedef  Proto::Box                              pr_box;
  typedef  Proto::DisjointBoxLayout                pr_dbl;
  template <int ncomp>
  class DataUtilities
  {
  public:
    typedef  Proto::LevelBoxData<double, ncomp>      pr_level_data;
    typedef  Proto::BoxData<     double, ncomp>      pr_box_data;


    ///
    /**
     **/
    shared_ptr<pr_dbl> getMetaData(const vector<pair<pr_box, int> > >&  a_layout)
    {
      shared_ptr<pr_dbl> retval;
      
      pr_mayday::Error("Doh! getMetaData not implemented")
      return retval;
    }

    ///
    /**
     **/
    void
    copyToHost( vector< double * >  &    a_host_data,
                const pr_level_data &    a_devi_data)
    {
      pr_mayday::Error("Doh! copyToHost not implemented")
    }

    ///
    /**
     **/
    void
    copyToDevice( pr_level_data             &    a_devi_data,
                  const vector< double * >  &    a_host_data)
    {
      pr_mayday::Error("Doh! copyToDevice not implemented")
    }

    ///
    /**
     **/
    void
    domainGhostReflectEven(pr_level_data             &    a_devi_data,
                           const pr_box              &    a_domain)
    {
      pr_mayday::Error("Doh!  ghostReflectEven not implemented");
    }

    ///
    /**
     **/
    void
    domainGhostReflectOdd(pr_level_data             &    a_devi_data,
                          const pr_box              &    a_domain)
    {
      pr_mayday::Error("Doh!  ghostReflectOdd not implemented");
    }

    ///
    /**
     **/
    template <int interpolation_order>
    void
    fillGhostOverCoarseFineBoundaries(pr_level_data             &    a_fine_data,
                                      const pr_level_data       &    a_coar_data,
                                      const pr_box              &    a_domain)
    {
      pr_mayday::Error("Doh!  fillGhostOverCoarseFineBoundaries not implemented");
    }
    
  };
}
#endif
