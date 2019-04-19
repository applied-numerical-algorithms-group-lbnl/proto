#include "EBProto.H"

using std::vector;

namespace Proto
{

/// Volume of Fluid Index
/**
*/
struct VolIndex
{
  inline
  VolIndex(const IntVect& a_ix,const  int& a_vofID)
  {
    define(a_ix, a_vofID);
  }


  inline void
  define(const IntVect& a_ix,const  int& a_vofID)
  {
    m_iv = a_ix;
    m_cellIndex = a_vofID;
  }
  
  inline
  bool operator<(const VolIndex& rhs) const
  {
    bool retval = false;
    if ( (m_iv < rhs.m_iv) || ((m_iv == rhs.m_iv) && (m_cellIndex < rhs.m_cellIndex))
    {
      retval = true;
    }
      return retval;
  }

  ///for stencil gymnastics
  void shift(const IntVect& a_iv)
  {
    m_iv += a_iv;
  }

  //cell of the VolIndex
  IntVect m_iv;

  //which VoF in the cell (0 to nVoF-1)
  int m_cellIndex;

};


///
/**
   FaceIndex is a very lightweight object used
   to distinguish faces.  It has two VolIndex (s)
   and a direction and a way to access them.
*/
struct FaceIndex
{
  ///
  /**
     if both cellindexes of the vofs are >= 0,
     then the face is interior.  otherwise the face is
     a boundary face.
   */
  FaceIndex(const VolIndex& a_vof1,
            const VolIndex& a_vof2,
            const int& a_direction)
  {
    define(a_vof1, a_vof2, a_direction);
  }
  FaceIndex(const VolIndex& a_vof1,
            const VolIndex& a_vof2,
            const int& a_direction)
  {
    if(a_vof1 < a_vof2)
    {
      m_vofLo = a_vof1;
      m_vofHi = a_vof2;
    }
    else
    {
      m_vofHi = a_vof1;
      m_vofLo = a_vof2;
    }
    m_direction = a_direction;
    m_isBoundary = ((m_vofLo.m_cellIndex < 0) ||(m_vofHi.m_cellIndex < 0));
  }



  ///
  /**
   */
  inline const bool& isBoundary() const
  {
    return  m_isBoundary;
  }

  ///for stencil gymnastics
  inline void shift(const IntVect& a_iv)
  {
    m_vofLo.shift( a_iv);
    m_vofHi.shift( a_iv);
  }

  VolIndex m_vofLo;
  VolIndex m_vofHi;
  int m_direction;
  bool m_isBoundary;
};


///
/**
   This is a list showing the connectivity of a given cell.  This
   has the special property of also being able to set itself to regular
   or covered by setting the current nodeimplem pointer to 1 or 0.
 */
class GraphCell
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
  
  GraphCell()
  {
    m_cellList = 0;
  }

  ~GraphCell()
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

  ///
  /**
     Get the faces in the direction and side
     for the vof in the list.  if the vof's cell
     index is not found in the list, abort.
     Use the input vof's grid index for a the grid index
     of the list.
   */
  inline vector<FaceIndex>
  getFaces(const VolIndex&       a_vof,
           const int&            a_idir,
           const Side::LoHiSide& a_sd,
           const ProblemDomain&  a_domain) const;

  ///
  inline vector<FaceIndex>
  getFaces(const IntVect&        a_iv,
           const int&            a_idir,
           const Side::LoHiSide& a_sd,
           const ProblemDomain&  a_domain) const;

  ///
  /**
     Return all the vofs in the list, using the input
     intvect for the gridIndex
   */
  inline vector<VolIndex>
  getVoFs(const IntVect& a_iv) const;

  ///
  /**
     Get all sets of connected vofs within the box (used in coarsening)
  */
  inline vector<vector<VolIndex> >  getVoFSets(const Box& a_box) const;

  ///
  /**
     The connectivity data at this point.
     If m_cellList == 0, node is covered.
     If m_cellList == 1, node is regular and parent is single-valued.
     otherwise, it is a real list (node may still be regular).
  */
  vector<IrrGraphNode>* m_cellList;

  private:
  GraphCell(const GraphNode& ebiin);
  void operator=(const GraphCell& ebiin);


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
    vector<VolIndex> getVoFs(const Point& a_iv) const;


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
    bool isRegular(const Box& a_box) const;

    ///
    bool isRegular(const Point& a_iv) const;

    ///
    bool isCovered(const Point& a_iv) const;


    ///
    bool isConnected(const VolIndex& a_vof1,
                     const VolIndex& a_vof2) const;

    ///
    /**
       Return the faces on the side and direction
       of the input VoF.
    */
    vector<FaceIndex> getFaces(const VolIndex& a_vof,
                               const int& a_idir,
                               const Side::LoHiSide& a_sd) const;


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


    inline bool isDefined() const
    {
      return m_isDefined;
    }
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
