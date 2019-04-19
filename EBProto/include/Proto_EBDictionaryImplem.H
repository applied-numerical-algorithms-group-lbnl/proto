
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
  EBDictionary(const GeometryService<order> & a_geoserv,
               const RealVect               & a_origin,
               const double                 & a_dx,
               const DisjointBoxLayout      & a_finestGrids,
               const Point                  & a_ghost,
               int  a_maxCoarsen)
  {
    shared_ptr<BaseIF> impfunc = a_geoserv.getBaseIF();
    m_grids.push_back(  a_finestGrids);
    m_domains.push_back(a_finestGrids.getDomain());

    int ilev = 0;
    while( ilev < a_maxCoarsen ) 
    {
      if(m_grids[ilev].coarsenable(2))
      {
        DisjointBoxLayout coarseGrids = grids[ilev];
        coarseGrids.coarsen(2);
        m_grids.push_back(  coarseGrids);
        m_domains.push_back(coarseGrids.getDomain());
        ilev++;
      }
      else
      {
        ilev = a_maxCoarsen+1;
      }
    }

    int numLev = m_grids.size();
    PR_assert(numLev >= 0);
    //resize the other vectors
    m_graphs.resize(numLev);
    m_voludata.resize(numLev);
    m_ebfadata.resize(numLev);
    m_registeredStencils.resize(numLev);
    for(int ifac = 0; ifac < DIM; ifac++) m_facedata[ifac].resize(numLev);

    //finest level of moment data comes from the geometry service
    //typically, you need graphs to be a bit bigger than the actual stencil
    m_geomgrow = a_ghost + Point::Ones(1);

    defineFinestLevel(a_geoserv, a_origin, a_dx);

    //coarser data is created via graph coarsening
    for(int lev = 1; ilev < numLev; ilev++)
    {
      PR_error("EBDictionary coarsening not yet implemented");
    }
    
  }    
    
  template<int order>
  void EBDictionary<order>::
  defineFinestLevel(const GeometryService<order> & a_geoserv,
                    const RealVect               & a_origin,
                    const double                 & a_dx )
  {
    template<int order>
    m_graphs[0].define(m_grids[0],   m_geomgrow);
    m_voludata[0].define(m_grids[0], m_geomgrow);
    m_ebfadata[0].define(m_grids[0], m_geomgrow);
    for(int ifac = 0; ifac < DIM; ifac++) m_facedata[ifac][0].define(m_grids[0], m_geomgrow);

    for(int ibox = 0; ibox < m_grids[0].size(); ibox++)
    {
      Box validBox = m_grids[0][ibox];
      Box grownBox = validBox.grow(geomgrow);
      Box interBox = grownBox & domains[0];
      if(impfunc->entireBoxRegular(interBox, a_origin, a_dx))
      {
        m_graphs[0][ibox].defineAsAllRegular(interBox, m_domains[0]);
      }
      else if(impfunc->entireBoxCovered(interBox, a_origin, a_dx))
      {
        m_graphs[0][ibox].defineAsAllCovered(interBox, m_domains[0]);
      }
      else
      {
        HostBoxData<int> regIrregCovered;
        vector< IrregNode<order> > nodes;
        a_geoserv.fillGRaph(regIrregCovered, nodes, validBox, interbox, m_domains[0]);
        vector<IrregGraph> irrgraph;
        getIrregGraphFromNodes<order>(irrgraph, nodes);

        m_graphs[0][ibox].defineFromPrimitives(regIrregCovered, irrgraph, intervox, m_domains[0]);
        //now fill in the moment information
        shared_ptr< vector<VolIndex>  > irregVoFs;
        m_graphs[0][ibox].getIrregVoFs( irregVoFs);
        m_voludata[0][ibox].define(     irregVoFs);
        for(int idir = 0; idir < DIM; idir++)
        {
          shared_ptr< vector<FaceIndex> > irregFaces;
          m_graphs[0][ibox].getIrregFaces( irregFaces, idir);
          m_facedata[idir][0][ibox].define(irregFaces);
        }
        for(int inode = 0; inode < nodes.size(); inode++)
        {
          const IrregNode<order>& node = nodes[inode];
          VolIndex vof(node.m_cell, node.m_cellIndex);
          m_voludata[0][ibox](vof, 0)  = node.m_volumeMoments;
          m_ebfadata[0][ibox](vof, 0)  = node.m_EBMoments;
          for(int idir = 0; idir < DIM; idir++)
          {
            for(SideIterator sit; sit.ok(); ++sit)
            {
              int facindex = ebp_index(idir, sit());
              vector<FaceIndex> faces = m_graphs[0][ibox].getFaces(vof, idir, sit());
              //finest level, single valued, etc
              PR_assert(faces.size() == 1);
              m_facdata[idir][0][ibox](faces[0], 0)  = node.m_faceMoments[facindex];
            }
          }
        }
      }
    }
  }

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
      m_registeredStencils[ilev][key] = shared_ptr<EBStencil>(new EBStencl(a_stencilName, a_domainBCName, a_ebbcName));
    }
  }


  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  template<int order>
  void EBDictionary<order>::
  close()
  {
    m_voludata = vector<shared_ptr<LevelData<VoluData> > >();
    m_ebfadata = vector<shared_ptr<LevelData<VoluData> > >();
    for(int idir = 0; idir < DIM; idir++)
      m_facedata[idir] = vector<shared_ptr<LevelData<FaceData> > >();
  }

  template<int order>
  shared_ptr<EBStencil>
  EBDictionary<order>::
  getEBStencil(const string & a_stencilName,
               const string & a_ebbcName, 
               const Box    & a_domain)
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
  }
  stencilKey key;
  key.m_ebbcName    = a_ebbcName;
  key.m_stencilName = a_stencilName;
  shared_ptr<EBStencil>   retval;
  if(m_registeredStencils[ilev].find(key) == m_registeredStencils[ilev].end())
  {
    PR_error("stencil not found");
  }
  else
  {
    retval = m_registeredStencils[ilev][key];
  }
}
#endif
