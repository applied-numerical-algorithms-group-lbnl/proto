#include "EBProto.H"

using std::vector;

namespace Proto
{
  template<CENTERING cent>
  class EBIndex
  {
  public:

    Point m_pt;
    //for faces, the high cell. for cells, this is the cell
    int   m_vofIDMe; 
    //for faces, the low cell. for cells, undefined
    int   m_vofIDLo;

    ///
    bool m_isBoundary;

    ///for vof construction
    EBIndex(Point a_pt, int a_cellInd) 
    {
      m_pt = a_pt;
      m_vofIDMe = a_cellInd;
    }

    //for sussing out a face's direction
    inline int direction() const
    {
      int retval = -1;
      if(cent == XFACE)
      {
        retval = 0;
      }
      else if(cent == YFACE)
      {
        retval = 1;
      }
      else if(cent == ZFACE)
      {
        retval = 2;
      }
      return retval;
    }
    //for getting the vofs of a face
    inline EBIndex<CELL> getVoF(Side::LoHiSide& a_sd)
    {
      EBIndex<CELL> retval;
      if(a_sd == Side::Lo)
      {
        int dir = direction();
        retval.m_pt = m_pt - Point::Basis(dir);
        retval.m_vofIDMe = m_vofIDLo;
      }
      else
      {
        retval.m_pt = m_pt;
        retval.m_vofIDMe = m_vofIDMe;
      }
      return retval;
    }

  };



  ///
  /**
     This is a list showing the connectivity of a given cell.  This
     has the special property of also being able to set itself to regular
     or covered by setting the current nodeimplem pointer to 1 or 0.
  */
  class CellNode
  {
   public:

    //each cell  contains a pointer to vector of these
    //the pointer has a special value for regular and covered cells.
    struct IrrGraphNode
    {
      //cell index is implied by place in the vector<IrrGraphNode>
      //grid index is implied by place in the HostData
      ///
      vector<int> m_arc[2*DIM];

    };
  
    void
    getFaceInfo(int   & a_vecsize,
                Point & a_mypt,
                int   & a_myVoFID,
                int   & a_otherVoFID
                bool  & a_isBoundary,
                const int            & a_dir,
                const  EBIndex<CELL>        & a_vof,
                const Side::LoHiSide & a_sd,
                const Box            & a_domain) const;

    CellNode()
    {
      m_cellList = 0;
    }

    ~CellNode()
    {
      if(isIrregular())
      {
        delete m_cellList;
      }
    }
    ///
    /**
       Return true if the node is covered (m_cellList==0).
    */
    inline bool isCovered() const
    {
      return m_cellList==0;
    }

    ///
    /**
       Return true if the node is regular.
    */
    inline bool isRegular() const;
    {
      return m_cellList==1;
    }

    ///
    /**
       Return true if the node is neither regular or covered.
    */
    inline bool isIrregular() const
    {
      return ((!isRegular()) && (!isCovered));
    }

    ///
    /**
       Set Node to regular.  If previously set to irregular,
       deletes their memory.
    */
    inline void defineAsRegular()
    {
      if(isIrregular()) delete m_cellList;
      m_cellList = 1;
    }


    //resets the size of the irregular list to 0
    inline void defineAsIrregular()
    {
      if(!isIrregular())
      {
        m_cellList = new vector<IrrGraphNode>();
      }
      m_cellList->resize(0);
    }

    ///
    /**
       Set Node to covered.  If previously set to irregular,
       deletes their memory.
    */
    inline void defineAsCovered()
    {
      if(isIrregular()) delete m_cellList;
      m_cellList = 1;
    }

    ///
    /**
       Get the faces in the direction and side
       for the vof in the list.  if the vof's cell
       index is not found in the list, abort.
       Use the input vof's grid index for a the grid index
       of the list.
    */
    inline vector<EBIndex<XFACE> >
    getXFaces(const  EBIndex<CELL> &       a_vof,
              const Side::LoHiSide& a_sd,
              const Box           & a_domain);


    inline vector<EBIndex<YFACE> >
    getYFaces(const  EBIndex<CELL> &       a_vof,
              const Side::LoHiSide& a_sd,
              const Box           & a_domain);


#if DIM==3
    inline vector<EBIndex<ZFACE> >
    getZFaces(const  EBIndex<CELL> &       a_vof,
              const Side::LoHiSide& a_sd,
              const Box           & a_domain);

#endif

    ///
    /**
       Return all the vofs in the list, using the input
       intvect for the gridIndex
    */
    inline vector< EBIndex<CELL> >
    getVoFs(const Point & a_iv) const;



    ///return vofs labelled as irregular (only places where an EB can happen)
    inline vector<EBIndex<BOUNDARY> >
    getCutFaces(const Point & a_iv) const;



    ///
    /**
       The connectivity data at this point.
       If m_cellList == 0, node is covered.
       If m_cellList == 1, node is regular and parent is single-valued.
       otherwise, it is a real list (node may still be regular).
    */
    vector<IrrGraphNode>* m_cellList;

  private:
    CellNode(const CellNode& ebiin);
    void operator=(const CellNode& ebiin);


  };


/// Geometric description within a box
/**
   EBGraph represents the geometric information
   of the domain at a given refinement
   within the boundaries of a particular box.
*/
  class EBGraph
  {
  public:
    ///
    /**
       The EBGraph can be all regular or all
       covered or contain irregular cells.
    */
    enum TAG
    {
      AllRegular = 0,
      AllCovered,
      HasIrregular
    };

    ///
    /**
       Makes an EBGraph whose contents are undefined.
    */
    EBGraph()
    {
      m_isDefined = false;
    }


    ///needed for leveldata
    EBGraph(const Box& a_grid)
    {
      m_isDefined = false;
    }

    ///
    ~EBGraph()
    {
    }

    ///
    /**
       Return the VoFs in the cell.
    */
    vector< EBIndex<CELL> > getVoFs(const Point& a_iv) const;


    ///
    /**
       Return true if every cell in the EBGraph is a regular cell.
    */
    bool isAllRegular() const
    {
      return (m_tag == AllRegular);
    }

    ///
    /**
       Return true neither all regular nor all covered
    */
    bool hasIrregular() const;

    ///
    /**
       Return true if  every cell in the EBGraph is a covered cell.
    */
    bool isAllCovered() const
    {
      return (m_tag == AllCovered);
    }
    ///
    /**
       Return true if a_iv is an irregular cell.
    */
    bool isIrregular(const Point& a_iv) const;

    ///
    /**
       Return true if every cell in a_box is a
       covered cell.
    */
    bool isCovered(const Box& a_box) const;

    ///
    /**
       Return true if every cell in a_box is a
       regular cell.
    */
    bool isRegular(const Box  & a_box) const;

    ///
    bool isRegular(const Point& a_iv) const;

    ///
    bool isCovered(const Point& a_iv) const;


    ///
    bool isConnected(const  EBIndex<CELL> & a_vof1,
                     const  EBIndex<CELL> & a_vof2) const;

    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    vector< EBIndex<XFACE> > getXFaces(const  EBIndex<CELL>      & a_vof,
                                       const Side::LoHiSide      & a_sd) const;


    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    vector< EBIndex<YFACE> > getYFaces(const  EBIndex<CELL>       & a_vof,
                                       const Side::LoHiSide      & a_sd) const;


#if DIM==3
    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    vector< EBIndex<ZFACE> > getZFaces(const  EBIndex<CELL>    & a_vof,
                                       const Side::LoHiSide    & a_sd) const;

#endif


    ///
    /**
       Completely construct the graph from IrregGraphs.
       This makes a DENSE representation of the graph.  To make
       sparse representations of all regular or all covered graphs,
       use setToAllRegular or setToAllCovered.\\
       irregGraph contains a complete list of irregular vofs. \\
       validRegion is the the region of the graph \\
       domain is the domain of compuation \\
       regIrregCovered = 1 for regular, 0 for irregular, -1 for covered.
       regIrregCovered must contain the valid region grown by 1. \\
       In chombospeak, \\
       regIrregCovered.box().contains(grow(validRegion, 1) & domain); \\
       otherwise an error will result.
       SingleValuedOnly signals that the geometry generator can only generate
       single-valued cut cells (GeometryService has this quality).  
       Coarsening can still create multi-valued cells.
    */
    void
    defineFromPrimitives(const HostBoxData<int>  &     a_regIrregCovered,
                         const vector<IrregGraph>&     a_irregGraph,
                         const Box               &     a_validRegion,
                         const Box               &     a_domain,
                         bool singleValuedOnly = true);

    ///
    /**
       Set the graph to all regular cells (sparse);
    */
    void defineAsAllRegular(const Box& a_region, const Box& a_domain);


    ///
    /**
       Set the graph to all regular cells (sparse);
    */
    void defineAsAllCovered(const Box& a_region, const Box& a_domain);

    ///
    const Box& getRegion() const
    {
      return m_region;
    }

    ///
    inline const Box& getDomain() const
    {
      return m_domain;
    }

    ///
    inline bool isDefined() const
    {
      return m_isDefined;
    }

    ////get the cut cells
    inline PointSet getIrregCells(const Box& a_box)
    {
      PR_assert(m_isDefined);

      PointSet irreg;
      if(hasIrregular())
      {
        irreg = m_internals->m_irregIVS;
        irreg &= a_box;
      }
      return irreg;
    }

    //get the cut vofs
    inline vector< EBIndex<CELL> > getIrregVoFs(const Box& a_box)
    {
      Box dombox = a_box & m_domain;
      PointSet irreg = getIrregCells(dombox);
      vector< EBIndex<Cell> > vofs;
      for(PointSetIterator ivsit(m_irreg);  ivsit.ok(); ++ivsit)
      {
        vector<EBIndex<CELL> > vofpt = m_graph.getVoFs(ivsit());
        vofs.insert(vofs.end(), vofpt.begin(), vofpt.end());
      }
      return vofs;
    }


    //get the cut faces (basically the same as above
    inline vector< EBIndex<BOUNDARY> > getCutFaces(const Box& a_box)
    {
      Box dombox = a_box & m_domain;
      PointSet irreg = getIrregCells(dombox);
      vector< EBIndex<BOUNDARY> > vofs;
      for(PointSetIterator ivsit(m_irreg);  ivsit.ok(); ++ivsit)
     {
        vector<EBIndex<CELL> > vofpt = m_graph.getCutFaces(ivsit());
        vofs.insert(vofs.end(), vofpt.begin(), vofpt.end());
      }
      return vofs;
    }

    ///get the faces where at least one vof is cut
    inline vector<EBIndex<XFACE> > getIrregXFaces(const Box& a_box)
    {
      //need to do this is a set at first to avoid double counting faces.
      std::unordered_set< EBIndex<XFACE> > resultSet;  
      Box dombox = a_box & m_domain;
      PointSet irreg = getIrregCells(dombox);
      for(PointSetIterator ivsit(m_irreg);  ivsit.ok(); ++ivsit)
      {
        vector<IrregCells> vofpt = m_graph.getVoFs(ivsit());
        for(int ivof = 0; ivof < vofpt.size(); ivof++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            vector<EBIndex< YFACE> > faces = getXFaces(vofpt[ivof],sit());
            for(int ifac = 0; ifac < faces.size(); ifac++)
            {
              resultSet.insert(faces[ifac]);
            }
          }
        }
      }
      
      vector< EBIndex<XFACE> > faces(resultSet.size());
      for (auto i = resultSet.begin(); i != resultSet.end();   ++i )
      {
        faces.push_back( *i );
      }
      
      return faces;
    }


    //get the faces where at least one vof is cut
    inline vector<EBIndex<YFACE> > getIrregYFaces(const Box& a_box)
    {
      //need to do this is a set at first to avoid double counting faces.
      std::unordered_set< EBIndex<YFACE> > resultSet;  
      PointSet irreg = getIrregCells(a_box);
      for(PointSetIterator ivsit(m_irreg);  ivsit.ok(); ++ivsit)
      {
        vector<EBIndex<CELL> > vofpt = m_graph.getVoFs(ivsit());
        for(int ivof = 0; ivof < vofpt.size(); ivof++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            vector<EBIndex< YFACE> > faces = getYFaces(vofpt[ivof],sit());
            for(int ifac = 0; ifac < faces.size(); ifac++)
            {
              resultSet.insert(faces[ifac]);
            }
          }
        }
      }
      
      vector< EBIndex<YFACE> > faces(resultSet.size());
      for (auto i = resultSet.begin(); i != resultSet.end();   ++i )
      {
        faces.push_back( *i );
      }
      
      return faces;
    }

#if DIM==3
    //get the faces where at least one vof is cut
    inline vector<EBIndex<ZFACE> > getIrregZFaces(const Box& a_box)
    {
      //need to do this is a set at first to avoid double counting faces.
      std::unordered_set< EBIndex<ZFACE> > resultSet;  
      PointSet irreg = getIrregCells(a_box);
      for(PointSetIterator ivsit(m_irreg);  ivsit.ok(); ++ivsit)
      {
        vector<EBIndex<CELL> > vofpt = m_graph.getVoFs(ivsit());
        for(int ivof = 0; ivof < vofpt.size(); ivof++)
        {
          for(SideIterator sit; sit.ok(); ++sit)
          {
            vector<EBIndex< YFACE> > faces = getZFaces(vofpt[ivof],sit());
            for(int ifac = 0; ifac < faces.size(); ifac++)
            {
              resultSet.insert(faces[ifac]);
            }
          }
        }
      }
      
      vector< EBIndex<ZFACE> > faces(resultSet.size());
      for (auto i = resultSet.begin(); i != resultSet.end();   ++i )
      {
        faces.push_back( *i );
      }
      
      return faces;
    }
#endif

  private:

    ///
    /**
       If this is allregular or allcovered,
       the HostBoxData below is undefined.
    */
    TAG m_tag;

    ///box over which this graph is defined
    Box m_region;

    ///computational domain at this level of refinement
    Box m_domain;

    ///
    bool m_isDefined;

    struct localData
    {

      ///
      /**
         the graph description of the box.
         If the pointer CellNode.m_cellList == 0, the cell is covered.
         If the pointer CellNode.m_cellList == 1, the cell is regular.
         Otherwise, the pointer goes to the node of the graph.
      */
      HostBoxData<CellNode, 1> m_graph;

      ///
      PointSet m_irregIVS;

    };


    // this is null in the case of m_tag == AllRegular || AllCovered
    shared_ptr<localData> m_internals;


  };
}

#include "Proto_EBGraphImplem.H"
#endif
