#ifndef _PROTO_LevelData_H_
#define _PROTO_LevelData_H_

#include "Proto_Point.H"
#include "Proto_Box.H"
#include "Proto_DisjointBoxLayout.H"
#include "Proto_WriteBoxData.H"
#include <cstdlib> //for size_t
#include <iostream>
#include <cstring> // for Writing data to .vtk files.

using std::shared_ptr;
namespace Proto 
{
///data over a disjointboxlayout with ghost cells.   
  template <class T>
  class LevelData
  {

  public:
    ///the int is the box index
    inline T & operator[](unsigned int a_index);

    ///the int is the box index
    inline const T & operator[](unsigned int a_index) const;

    ///
    inline LevelData();

    /// 
    inline LevelData(  const DisjointBoxLayout& a_grids, const Point& a_grow);

    /// 
    inline void define(const DisjointBoxLayout& a_grids, const Point& a_grow);

    ///
    inline unsigned int size() const;

    ///
    inline void setToZero();

    /// fills ghost cells with valid data
    inline void exchange();
  
    ///fills ONLY VALID DATA
    /** call exchange if you want ghost cells filled */
    inline void copyTo(LevelData<T>& a_dest) const;

    /// gets the dbl PC: added.
    inline DisjointBoxLayout getDBL(){return m_grids;};

  private: 

    vector<shared_ptr<T> >        m_data;
    Point                         m_grow;
    DisjointBoxLayout             m_grids;
    bool                          m_isDefined;
  };
  /// Stand-alone function that writes a LevelData on a domain-decomposed rectangle into a .vtk file.
  inline 
  void WriteData(LevelData<BoxData<double> >& a_phi,int a_iter,double a_dx,std::string a_str);
} //end namespace Proto
#include "Proto_LevelDataI.H"
#endif
