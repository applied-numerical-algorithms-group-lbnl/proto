
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
  EBDictionary::
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

    int numLev = m_grids.size();
    PR_assert(numLev >= 0);
    //resize the other vectors
    m_voludata.resize(numLev);
    m_registeredStencils.resize(numLev);
    for(int ifac = 0; ifac < 2*DIM; ifac++) m_facedata[ifac].resize(numLev);

    //finest level of moment data comes from the geometry service
    //typically, you need graphs to be a bit bigger than the actual stencil
    Point geomgrow = a_ghost + Point::Ones(1);

    m_graphs[0].define(m_grids[0], geomgrow);
    m_voludata[0].define(m_grids[0], geomgrow);
    for(int ifac = 0; ifac < 2*DIM; ifac++) m_facedata[ifac][0].define(m_grids[0], geomgrow);



    //coarser data is created via graph coarsening
    for(int lev = 1; ilev < numLev; ilev++)
    {
      PR_error("EBDictionary coarsening not yet implemented");
    }
    
    
  }


  inline void registerStencil(string a_stencilName)
  {
  }


  ///after this is done, you cannot register any more stencils.  This will delete the moment information.
  inline void close();

  inline shared_ptr<EBStencil> getEBStencil(const string & a_stencilName, 
                                            const Box    & a_domain);
protected:

  ///these vectors are of the length of number of coarsenings
  ///
  vector<shared_ptr<vector<VoluData> > >       m_voludata;

  vector<Box> m_domains; //keys to getting to what level is being called
  ///
  vector<shared_ptr<vector<FaceData> > >       m_facedata[2*DIM];

  ///
  vector< map<shared_ptr<EBStencil>, string> >   m_registeredStencils;
  
  vector<DisjointBoxLayout> m_grids;
private:
  EBDictionary();
  EBDictionary(  const EBDictionary<order>& a_inp);
  void operator=(const EBDictionary<order>& a_inp);


};
}
#endif
