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
    typedef  Proto::BoxData<double, 1>               pr_boxdata_sca;
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
    typedef  Proto::BoxData<     double, ncomp>      pr_box_data_ncomp;
    typedef  Proto::BoxData<     double,  1   >      pr_box_data_sca;
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
      int a_domain_face_dir,  int    a_zero_lo_one_hi,  int a_sten_sign,
      int a_force_num_ghost = -1)  
    {
      auto dit   = a_devi_data.dataIterator();
      auto grids = a_devi_data.layout();
      auto dombx = grids.domain().box();
      bool initToZero = true;
      int num_ghost = a_devi_data.ghost()[a_domain_face_dir];
      if(a_force_num_ghost > 0)
      {
        num_ghost = a_force_num_ghost;
      }
      
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        
        pr_box              valid = grids[dit[ibox]];
        vector<pr_sten>  stencils; 
        vector<pr_box>   appboxes;

        getDomainBCReflectionStencils(stencils, appboxes, valid, dombx, num_ghost,
                                      a_domain_face_dir, a_zero_lo_one_hi, a_sten_sign);
        // these vectors are empty away from the domain boundary
        for(int isten = 0; isten < stencils.size(); isten++)
        {
          for(int icomp = 0; icomp < ncomp; icomp++)
          {
            pr_box_data_ncomp& allphi  = a_devi_data[dit[ibox]];
            pr_box_data_sca    phicomp = Proto::slice(allphi, icomp);
            stencils[isten].apply(phicomp, phicomp, appboxes[isten], initToZero, 1.0);
          }
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
      int a_domain_face_dir,  int    a_zero_lo_one_hi,  
      int a_force_num_ghost = -1)  
    {
      int sten_sign = 1;
      domainGhostReflection(a_devi_data,
                            a_domain, a_domain_face_dir,
                            a_zero_lo_one_hi, sten_sign,
                            a_force_num_ghost);
    }

    ///
    /**
       fill ghost cells outside domain with an odd reflection of interior data
     **/
    static void
    domainGhostReflectOdd(
            pr_level_data             &    a_devi_data,
            const pr_box              &    a_domain,
            int a_domain_face_dir,  int a_zero_lo_one_hi,
            int a_force_num_ghost = -1)  
    {
      int sten_sign = -1;
      domainGhostReflection(a_devi_data,
                            a_domain, a_domain_face_dir,
                            a_zero_lo_one_hi, sten_sign,
                            a_force_num_ghost);
    }

    /// also calls exchange
    static void
    enforceReflectionBCs(pr_level_data &    a_devi_data,
                         const pr_box  &    a_domain,
                         int                a_reflect_sign_lo[DIM],
                         int                a_reflect_sign_hi[DIM],
                         int a_force_num_ghost = -1)  
    {
      const int zero_lo = 0;  const int one_hi = 1; 
      for(int idir = 0; idir < DIM; idir++)
      {
        domainGhostReflection(a_devi_data, a_domain, idir, zero_lo, a_reflect_sign_lo[idir], a_force_num_ghost);
        domainGhostReflection(a_devi_data, a_domain, idir,  one_hi, a_reflect_sign_hi[idir], a_force_num_ghost);
      }
    }

    ///
    /**
     *  \param a_fine_data - data whose ghost is getting filled
     *  \param a_coar_data - data   at the next coarser level
     *  \param a_domain_coar    - domain at the next coarser level
     *  \param ref_ratio   - refrat to the next coarser level
     *  \param reflect_sign_lo - reflection sign for the low side of the domain.
     *  \param reflect_sign_hi - reflection sign for the high side of the domain.

     
    Important caveats:
     *  --- For now, I am  only dealing with second order operators.
     *  --- I need to  fill only one ghost cell over the coarse-fine interface.  
     *  --- To top it all off, I can assume that the grids are properly nested.
     *  --- These special cases make this process a lot easier.
     *  --- The BC stuff is in here because it is an easy way for me to get
     *      off the ground and know that interpolation stencil will always
     *      have valid data.  This is probably a bad idea for the long term
     *      as we may have convergence problems for interpolation order > 2.
     * --- A lower dimension extrapolation stencil at each domain face is
     *      probably a better idea and (with some work) could be automated.
     * --- I am ignoring  interpolation_order for now since I will just be using Quadratic.
     *     Higher order operators should do something higher order at coarse-fine interfaces.
     * --- Implementing higher order coarse-fine interpolation requires that the data 
     *     has at least enough ghost cells to satisfy the span of the stencil.
     *     The only way I can think of unmaking this ghost-cell/order connection is 
     *     by using an awful lot of communications (repeated exchanges).   
     * --- Yes, the ghost cell computation is complicated (but in a good way) by the refinement ratio.
     *     (you are guaranteed nref ghost cells because you are guaranteed one coarse cell).
     dtg 
     January, 2024
    **/
    template <int interpolation_order>
    static void
    fillGhostOverCoarseFineBoundaries(
            pr_level_data        & a_data_fine,
            const pr_level_data  & a_data_coar,
            const pr_box         & a_domain_coar,
            const pr_pt          & a_ref_ratio,
            int                    a_reflect_sign_lo[DIM],
            int                    a_reflect_sign_hi[DIM])  
    {
      ///always
      const bool overwrite = true; const double applyScale = 1.;
      
      /**
         Get coarse data on coarsened fine layout.
      **/
      typedef Proto::InterpStencil<double> pr_interp;
      pr_interp quad_sten = pr_interp::Quadratic(a_ref_ratio);
      
      pr_pt ghost       = a_data_fine.ghost();
      pr_dbl  dbl_fine  = a_data_fine.layout();
      pr_dbl  dbl_cofi  =    dbl_fine.coarsen(a_ref_ratio);
      
      pr_level_data  data_cofi(dbl_cofi, ghost);
      pr_level_data  copy_fine(dbl_fine, pr_pt::Zeros());///just valid data
      a_data_fine.copyTo(copy_fine);
      a_data_coar.copyTo(data_cofi);
      //this also does exchange
      enforceReflectionBCs(data_cofi, a_domain_coar, a_reflect_sign_lo, a_reflect_sign_hi);

      auto dit   = a_data_fine.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto         & fab_fine = a_data_fine[dit[ibox]];
        const auto   & fab_cofi =   data_cofi[dit[ibox]];
        pr_box  valid_fine      =    dbl_fine[dit[ibox]];
        pr_box  grown_fine = valid_fine.grow(pr_pt::Ones());
        for( int icomp = 0; icomp < ncomp ; icomp++)
        {
          pr_box_data_sca sca_fine = Proto::slice(fab_fine, icomp);
          pr_box_data_sca sca_cofi = Proto::slice(fab_cofi, icomp);
          ///interp_stencil follows the chombo standard of dest, src so sca_fine is the output
          quad_sten.apply(sca_fine, sca_cofi, grown_fine, overwrite, applyScale);
        }
      }
      /**
         Copy back valid data.   This has no ghost cells so C/F cells should not be affected.
      **/
      copy_fine.copyTo(a_data_fine);
    }

    ///
    /**
       same as above where all the coarse data is zero
     **/
    template <int interpolation_order>
    static void
    homogeneousCFInterp(
            pr_level_data        & a_data_fine,
            const pr_pt          & a_ref_ratio)
    {
      ///always
      const bool overwrite = true; const double applyScale = 1.;
      
      /**
         Get coarse data on coarsened fine layout.
      **/
      typedef Proto::InterpStencil<double> pr_interp;
      pr_interp quad_sten = pr_interp::Quadratic(a_ref_ratio);
      
      pr_pt ghost       = a_data_fine.ghost();
      pr_dbl  dbl_fine  = a_data_fine.layout();
      pr_dbl  dbl_cofi  =    dbl_fine.coarsen(a_ref_ratio);
      
      pr_level_data  data_cofi(dbl_cofi, ghost);
      data_cofi.setVal(0.);
      pr_level_data  copy_fine(dbl_fine, pr_pt::Zeros());///just valid data
      auto dit   = a_data_fine.dataIterator();
      for(int ibox = 0; ibox < dit.localSize(); ibox++)
      {
        auto         & fab_fine = a_data_fine[dit[ibox]];
        const auto   & fab_cofi =   data_cofi[dit[ibox]];
        pr_box  valid_fine      =    dbl_fine[dit[ibox]];
        pr_box  grown_fine = valid_fine.grow(pr_pt::Ones());
        for( int icomp = 0; icomp < ncomp ; icomp++)
        {
          pr_box_data_sca sca_fine = Proto::slice(fab_fine, icomp);
          pr_box_data_sca sca_cofi = Proto::slice(fab_cofi, icomp);
          ///interp_stencil follows the chombo standard of dest, src so sca_fine is the output
          quad_sten.apply(sca_fine, sca_cofi, grown_fine, overwrite, applyScale);
        }
      }
      /**
         Copy back valid data.   This has no ghost cells so C/F cells should not be affected.
      **/
      copy_fine.copyTo(a_data_fine);
    }
    
  };
}
#endif
