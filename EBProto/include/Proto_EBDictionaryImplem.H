
#ifndef _Proto_EBDictionaryImplem_H_
#define _Proto_EBDictionaryImplem_H_



namespace Proto
{
  ///
  /**
     a_ghost is the number of ghost cells your data will hold.
     It must be at least as big as the stencil span in each direction..
  */
  template <int order>
  EBDictionary<order>::
  EBDictionary(const shared_ptr<GeometryService<order> >    & a_geoserv)
  {
    m_geoserv = a_geoserv;
  }    
    
  ///
  template<int order>
  void EBDictionary<order>::
  registerStencil(string a_stencilName,
                  string a_domainBCName[2*DIM],
                  string a_ebbcName)
  {
    stencilKey key;
    key.m_ebbcName    = a_ebbcName;
    key.m_stencilName = a_stencilName;
    for(int ilev = 0; ilev < m_grids.size(); ilev++)
    {
      m_registeredStencils[ilev][key] = shared_ptr<EBStencil>(new EBStencil(a_stencilName, a_domainBCName, a_ebbcName, m_grids[ilev], m_dataGrow));
    }
  }


  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  template<int order>
  void EBDictionary<order>::
  close()
  {
    m_isClosed = true;
    m_voludata = vector<shared_ptr<LevelData<VoluData> > >();
    m_ebfadata = vector<shared_ptr<LevelData<VoluData> > >();
    for(int idir = 0; idir < DIM; idir++)
    {
      m_facedata[idir] = vector<shared_ptr<LevelData<FaceData> > >();
    }
  }
  ///
  template<int order>
  shared_ptr<EBStencil>
  EBDictionary<order>::
  getEBStencil(const string & a_stencilName,
               const string & a_ebbcName, 
               const Box    & a_domain,
               const int    & a_boxid)
  {
    int level = -1;
    bool found = false;
    for(int ilev = 0; ilev < m_domains.size(); ilev++)
    {
      if(m_domains[ilev]== a_domain)
      {
        level = ilev;
        found = true;
        break;
      }
    }
    PROTO_ASSERT(found, "could not find domain in getEBStencil");
    stencilKey key;
    key.m_ebbcName    = a_ebbcName;
    key.m_stencilName = a_stencilName;
    shared_ptr<EBStencil>   retval;
    if(m_registeredStencils[ilev][a_boxid].find(key) == m_registeredStencils[ilev][a_boxid].end())
    {
      PR_error("stencil not found");
    }
    else
    {
      retval = m_registeredStencils[ilev][a_boxid][key];
    }
  }
}
#endif
