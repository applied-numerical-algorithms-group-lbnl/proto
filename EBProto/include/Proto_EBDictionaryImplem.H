
#ifndef _Proto_EBDictionaryImplem_H_
#define _Proto_EBDictionaryImplem_H_



namespace Proto
{
  ///
  /**
     a_ghost is the number of ghost cells your data will hold.
     It must be at least as big as the stencil span in each direction..
  */
  template <int order, typename T, CENTERING srcCenter, CENTERING dstCenter>
  EBDictionary<order, T, srcCenter, dstCenter>::
  EBDictionary(const shared_ptr<GeometryService<order> >    & a_geoserv,
               const DisjointBoxLayout    a_grids,
               Point  a_srcGhost,  Point  a_dstGhost)
  {
    m_geoserv  = a_geoserv;
    m_grids    = a_grids;
    m_srcGhost = a_srcGhost;
    m_dstGhost = a_dstGhost;
    m_isClosed = false;
    m_registeredStencils.resize(m_grids.size());
  }    
    
  ///
  template <int order, typename T, CENTERING srcCenter, CENTERING dstCenter>
  void
  EBDictionary<order, T, srcCenter, dstCenter>::
  registerStencil(string a_stencilName,
                  string a_domainBCName[2*DIM],
                  string a_ebbcName)
  {
    PR_assert(!m_isClosed);

    stencilKey key;
    key.m_ebbcName    = a_ebbcName;
    key.m_stencilName = a_stencilName;
    for(int ibox = 0; ibox < m_grids.size(); ibox++)
    {
HERE NEED TO SEND IN GEOMETRY INFO
      m_registeredStencils[ibox][key] = 
        shared_ptr<>(new ebstencil_t(a_stencilName, a_domainBCName, a_ebbcName, m_grids[ibox], m_src));
    }
  }


  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  template <int order, typename T, CENTERING srcCenter, CENTERING dstCenter>
  void
  EBDictionary<order, T, srcCenter, dstCenter>::
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
  template <int order, typename T, CENTERING srcCenter, CENTERING dstCenter>
  shared_ptr<EBStencil>
  EBDictionary<order, T, srcCenter, dstCenter>::
  getEBStencil(const string & a_stencilName,
               const string & a_ebbcName, 
               const int    & a_boxid)
  {
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
