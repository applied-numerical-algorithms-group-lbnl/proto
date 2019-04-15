#include "EBProto.H"

using std::vector;

namespace Proto
{
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
       Completely construct the graph from IrregNodes.
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
       otherwise an error will result
    */
    void
    defineFromPrimitives(const HostBoxData<int> &     a_regIrregCovered,
                         const vector<IrregNode>&     a_irregGraph,
                         const Box              &     a_validRegion,
                         const Box              &     a_domain);

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
         If the pointer GraphNode.m_thisNode == 0, the cell is covered.
         If the pointer GraphNode.m_thisNode == 1, the cell is regular.
         Otherwise, the pointer goes to the node of the graph.
      */
      HostBoxData<GraphNode> m_graph;

      ///
      PointSet* m_irregIVS;

      ///
      PointSet* m_multiIVS;
    };


    //this is to make this a ref-counted object
    shared_ptr<localData> m_internals;


  };
}

#include "Proto_EBGraphImplem.H"
#endif
