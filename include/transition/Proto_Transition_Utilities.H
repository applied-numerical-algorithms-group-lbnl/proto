#ifndef __Proto_Transition_Utilities_H__
#define __Proto_Transition_Utilities_H__

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
  ///class has no data
  class MetaUtilities
  {
  public:
    typedef  Proto::MayDay<void>                     pr_mayday;
    typedef  Proto::Point                            pr_pt;
    typedef  Proto::Box                              pr_box;
    typedef  Proto::DisjointBoxLayout                pr_dbl;
    typedef  Proto::ProblemDomain                    pr_dom;
    typedef  Proto::BoxPartition                     pr_partition;
    
    ///Get a Proto::DisjointBoxLayout that conforms to the input meta data
    /** 
        A box here is its point refined by the box size.   
        pair.first = the valid box of the layout.
        pair.second = the processor assignment of said box.
        The input layout can be in any order or have any distribution but 
        the boxes MUST BE DISJOINT.
    **/
    static shared_ptr< pr_dbl > getLayout(
      const vector<pair<pr_pt, unsigned int> >   &  a_layout,
      const pr_pt                                &  a_boxSize,
      const pr_dom                                  &  a_domain)
    {

      // Input Checks
      PROTO_ASSERT(a_layout.size() > 0,
                   "DataUtilities::getLayout | Error: empty input");
      PROTO_ASSERT(a_domain.box().coarsenable(a_boxSize),
                   "DataUtilities::getLayout | Error: domain is not coarsenable by the given box size.");
      pr_box patchDomain = a_domain.box().coarsen(a_boxSize);
      for (auto iter : a_layout)
      {
        pr_pt patch = iter.first;
        PROTO_ASSERT(patch.min() >= 0,
                     "DataUtilities::getLayout | Error: Patches must be in the non-negative index space");
        PROTO_ASSERT(patchDomain.contains(patch),
                     "DataUtilities::getLayout | Error: Patch detected outside of the given domain");
      }

      // Box Partition
      std::shared_ptr<pr_partition> partition = std::make_shared<pr_partition>(a_domain.coarsen(a_boxSize), a_layout);
      // Box Layout
      shared_ptr< pr_dbl>  retval(new pr_dbl());
      retval->define(partition, a_boxSize);
      return retval;
    }
  };
  
  ///This class also has no data. Just a pile of static functions, really.
  template <int ncomp>
  class DataUtilities
  {
  public:
    typedef  Proto::MayDay<void>                     pr_mayday;
    typedef  Proto::Point                            pr_pt;
    typedef  Proto::Box                              pr_box;
    typedef  Proto::DisjointBoxLayout                pr_dbl;
    typedef  Proto::ProblemDomain                    pr_dom;
    typedef  Proto::BoxPartition                     pr_partition;
    typedef  Proto::LevelBoxData<double, ncomp>      pr_level_data;
    typedef  Proto::BoxData<     double, ncomp>      pr_box_data;
    typedef  Proto::Stencil<     double>             pr_sten;


    ///
    static void
    copyToDevice(double       * a_devi_data,
                 const double * a_host_data,
                 size_t a_dat_size)
    {
#ifdef PROTO_ACCEL      
      Proto::proto_memcpy<Proto::HOST, Proto::DEVICE>(a_host_data, a_devi_data, a_dat_size);
#else      
      Proto::proto_memcpy<Proto::HOST, Proto::HOST  >(a_host_data, a_devi_data, a_dat_size);
#endif      
    }
    
    ///
    static void
    copyToHost(double       * a_host_data,
               const double * a_devi_data,
               size_t a_dat_size)
    {
#ifdef PROTO_ACCEL      
      Proto::proto_memcpy<Proto::DEVICE, Proto::HOST>(a_devi_data, a_host_data, a_dat_size);
#else      
      Proto::proto_memcpy<Proto::HOST  , Proto::HOST>(a_devi_data, a_host_data, a_dat_size);
#endif      
    }
    ///
    /**
     *  \param a_host_data - one buffer per box
     **/ 
    static void
    copyToHost(
            vector< double * >  &    a_host_data,
            const pr_level_data &    a_devi_data)
    {
      auto dit = a_devi_data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        const double* devi_ptr = a_devi_data[dit[ibox]].data();
        double      * host_ptr = a_host_data[ibox];
        const pr_box& fabbox   = a_devi_data[dit[ibox]].box();
        size_t num_points      = fabbox.size();
        size_t data_size       = ncomp*sizeof(double)*num_points;
        copyToHost(host_ptr, devi_ptr, data_size);
      }
      
    }


    ///
    /**
     *  \param a_devi_data  device data to be filled in
     *  \param a_host_data  matching host pointers for each BoxData in a_devi_data
     *  The host data has to be in fortran ordering and must contain the same footprint 
     *  as the device BoxData (it has to include ghost cell memory).
     **/
    static void
    copyToDevice(
            pr_level_data             &    a_devi_data,
            const vector< double * >  &    a_host_data) 
    {

      auto dit = a_devi_data.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        double*       devi_ptr = a_devi_data[dit[ibox]].data();
        const double* host_ptr = a_host_data[ibox];
        const pr_box& fabbox   = a_devi_data[dit[ibox]].box();
        size_t num_points      = fabbox.size();
        size_t data_size       = ncomp*sizeof(double)*num_points;
        copyToDevice(devi_ptr, host_ptr, data_size);
      }


    }


    static void
    getDomainBCReflectionStencils(vector<pr_sten>  & a_stencils,
                                  vector<pr_box>   & a_appboxes,
                                  const  pr_box    & a_validBox,
                                  const  pr_box    & a_domain,
                                  int                a_num_ghost,
                                  int                a_idir,
                                  int                a_zero_lo_one_hi,
                                  int a_sten_sign )//1 reflect even, -1 reflect odd
    
    {
      pr_pt   lobox = a_validBox.low();
      pr_pt   hibox = a_validBox.high();
      pr_pt   lodom = a_domain.low();
      pr_pt   hidom = a_domain.high();
      pr_box srcBox = a_validBox.grow(pr_pt::Ones(a_num_ghost));
      a_stencils.resize(   0);
      a_appboxes.resize(0);
      double  stenSign = double(a_sten_sign);

      if(a_zero_lo_one_hi == 0) //low side
      {
        pr_pt loapp  = srcBox.low();
        pr_pt hiapp  = srcBox.high();
        loapp[a_idir]  = lodom[a_idir]-1;
        hiapp[a_idir]  = lodom[a_idir]-1;

        if(lobox[a_idir] == lodom[a_idir])
        {
          a_stencils.resize(a_num_ghost);
          a_appboxes.resize(a_num_ghost);

          pr_box appBox0(loapp, hiapp);
          for(int igh = 0; igh < a_num_ghost; igh++)
          {
            a_appboxes[igh] = appBox0.shift(a_idir, -igh);
            a_stencils[igh] = (stenSign)*Proto::Shift(pr_pt::Basis(a_idir, 2*igh+1));
          }
        }
      }
      else if (a_zero_lo_one_hi == 1)// high side
      {
        if(hibox[a_idir] == hidom[a_idir]) 
        {
          a_stencils.resize(a_num_ghost);
          a_appboxes.resize(a_num_ghost);

          pr_pt loapp = srcBox.low();
          pr_pt hiapp = srcBox.high();
          loapp[a_idir] = hidom[a_idir]+1;
          hiapp[a_idir] = hidom[a_idir]+1;

          pr_box appBox0(loapp, hiapp);
        
          for(int igh = 0; igh < a_num_ghost; igh++)
          {
            a_appboxes[igh] = appBox0.shift(a_idir, igh);
            a_stencils[igh] = (stenSign)*Proto::Shift(pr_pt::Basis(a_idir, -2*igh-1));
          }
        }
      
      }
      else
      {
        pr_mayday::Error("Doh! getDomainBCs: bogus zero_lo_one_hi (must be 0 or 1)");
      }
    }


    ///
    /**
     **/
    static void
    domainGhostReflection(
      pr_level_data             &    a_devi_data,
      const pr_box              &    a_domain,
      int a_domain_face_dir,  int    a_zero_lo_one_hi,  int a_sten_sign)  
    {
      auto dit   = a_devi_data.dataIterator();
      auto grids = a_devi_data.layout();
      auto dombx = grids.domain().box();
      bool initToZero = true;
      int num_ghost = a_devi_data.ghost()[a_domain_face_dir];
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        
        pr_box              valid = grids[dit[ibox]];
        vector<pr_sten>  stencils; 
        vector<pr_box>   appboxes;

        getDomainBCReflectionStencils(stencils, appboxes, valid, dombx, num_ghost,
                                      a_domain_face_dir, a_zero_lo_one_hi, a_sten_sign);
        // these vectors are empty away from the domain boundary
        auto& datafab = a_devi_data[dit[ibox]];
        for(int isten = 0; isten < stencils.size(); isten++)
        {
          stencils[isten].apply(datafab, datafab, appboxes[isten], initToZero, 1.0);
        }
      }
    }

    ///
    /**
       fill ghost cells outside domain with an even reflection of interior data
     **/
    static void
    domainGhostReflectEven(
      pr_level_data             &    a_devi_data,
      const pr_box              &    a_domain,
      int a_domain_face_dir,  int    a_zero_lo_one_hi)  
    {
      int sten_sign = 1;
      domainGhostReflection(a_devi_data,
                            a_domain, a_domain_face_dir,
                            a_zero_lo_one_hi, sten_sign);
    }

    ///
    /**
       fill ghost cells outside domain with an odd reflection of interior data
     **/
    static void
    domainGhostReflectOdd(
            pr_level_data             &    a_devi_data,
            const pr_box              &    a_domain,
            int a_domain_face_dir,  int a_zero_lo_one_hi)
    {
      int sten_sign = -1;
      domainGhostReflection(a_devi_data,
                            a_domain, a_domain_face_dir,
                            a_zero_lo_one_hi, sten_sign);
    }

    ///
    /**
     *  \param a_fine_data - data whose ghost is getting filled
     *  \param a_coar_data - data   at the next coarser level
     *  \param a_domain    - domain at the next coarser level
     *  \param ref_ratio   - refrat to the next coarser level
     **/
    template <int interpolation_order>
    static void
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
