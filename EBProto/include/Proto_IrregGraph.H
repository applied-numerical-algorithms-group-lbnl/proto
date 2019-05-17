
#ifndef _Proto_IRREGGRAPH_H_
#define _Proto_IRREGGRAPH_H_


#include "Proto_Point.H"

using std::vector;
namespace Proto
{


  ///simple struct to allow me to keep EBGraph untemplated
  struct IrregGraph
  {
  public:

    ///
    Point     m_cell;
    int       m_cellIndex;

    ///
    /**
       Indicies into a_nodes to show connectivity.
       If the arc is to  an irregular
       cell, the index is the unique index of the vof
       in the cell.
       For arcs to regular cells, the arc = -2
       If the arc is to the domain boundary, arc  = -1.
    */
    vector<int>      m_arc[2*DIM];


  };

}


#endif
