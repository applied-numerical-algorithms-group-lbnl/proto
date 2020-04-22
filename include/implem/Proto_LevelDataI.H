#ifndef _PROTO_LevelDataImplem_H_
#define _PROTO_LevelDataImplem_H_

namespace Proto
{

  ///the int is the box index
  template<class T>
  T & 
  LevelData<T>::
  operator[](unsigned int a_index)
  {
    
    PROTO_ASSERT(m_isDefined,"trying to access undefined leveldata");
    PROTO_ASSERT(a_index < m_data.size(),"bogus index sent to leveldata");
    DisjointBoxLayout::activeBox(a_index);

    return (*(m_data[a_index]));
  }

  ///the int is the box index
  template<class T>
  const T & 
  LevelData<T>::
  operator[](unsigned int a_index) const
  {
    
    PROTO_ASSERT(m_isDefined,"trying to access undefined leveldata");
    PROTO_ASSERT(a_index < m_data.size(),"bogus index sent to leveldata");
    DisjointBoxLayout::activeBox(a_index);
    
    return (*(m_data[a_index]));
  }

  ///
  template<class T>
  LevelData<T>::
  LevelData()
  {
    m_isDefined = false;
  }


  ///
  template<class T>
  LevelData<T>::
  LevelData(  const DisjointBoxLayout& a_grids, const Point& a_grow)
  {
    define(a_grids, a_grow);
  }

  ///
  template<class T>
  void 
  LevelData<T>::
  define(const DisjointBoxLayout& a_grids, const Point& a_grow)
  {
    m_isDefined = true;
    m_grids = a_grids;
    m_grow  = a_grow;

    m_data.resize(m_grids.size());
    for(unsigned int ibox = 0; ibox < m_grids.size(); ibox++)
    {
      Box dblbox  = m_grids[ibox];
      Box databox = dblbox.grow(a_grow);
      
      T*  newdata = new T(databox, false);

      m_data[ibox] = shared_ptr<T>(newdata);
    }
  }

  ///
  template<class T>
  unsigned int 
  LevelData<T>::size() const
  {
    return m_data.size();
  }

  ///
  template<class T>
  void   
  LevelData<T>::setToZero()
  {
    for(unsigned int ibox = 0; ibox < m_data.size(); ibox++)
    {
      m_data[ibox]->setVal(0);
    }
  }

  /// fills ghost cells with valid data
  template<class T>
  void 
  LevelData<T>::exchange()
  {
    if(m_grow != Point::Zeros())
    {
      //loop through my valid data
      for(unsigned int sbox = 0; sbox < m_data.size(); sbox++)
      {
        Box validBoxSrc = m_grids[sbox];
        vector<NeighborElement> neighbors = m_grids.getNeighbors(sbox);
        for(unsigned int ineigh = 0; ineigh < neighbors.size(); ineigh++)
        {
          //neighbors includes the current box
          unsigned int dbox = neighbors[ineigh].m_boxID;
          Point pshift =      neighbors[ineigh].m_periodicShift;
          if((dbox  != sbox) || (pshift != Point::Zeros()))
          {
            Box grownNeighborBox = (*m_data[dbox]).box();

            Box shiftedBoxDst = grownNeighborBox.shift(-pshift);
            Box intersectSrc = validBoxSrc & shiftedBoxDst;
            if(!intersectSrc.empty())
            {
              m_data[sbox]->copyTo(*m_data[dbox], intersectSrc, pshift);
            }
          }
        }
      }
    }
  }
  
  ///fills ONLY VALID DATA
  /** call exchange if you want ghost cells filled */
  template<class T>
  void 
  LevelData<T>::
  copyTo(LevelData<T>& a_dest) const
  {
    if(a_dest.m_data != m_data)
    {
      if(a_dest.m_grids == m_grids)
      {
        //loop through my valid data
        for(unsigned int sbox = 0; sbox < m_data.size(); sbox++)
        {
          Box validBoxSrc = m_grids[sbox];
          m_data[sbox]->copyTo(*a_dest.m_data[sbox], validBoxSrc);
        }
      }
      else
      {
        //just implementing the n^2 algorithm here for now-- should probably do moving window
        //but this might never get called
        for(unsigned int sbox = 0; sbox < m_data.size(); sbox++)
        {
          Box validBox = m_grids[sbox];
          for(unsigned int dbox = 0; dbox < a_dest.m_data.size(); dbox++)
          {
            Box destbox = a_dest.m_grids[dbox];
            Box intersect = validBox & destbox;
            if(!intersect.empty())
            {
              m_data[sbox]->copyTo(*a_dest.m_data[dbox], intersect);
            }
          }
        }
      }
    }
  }
  void WriteData(LevelData<BoxData<double> >& a_phi,int a_iter,double a_dx,std::string a_str)
  {
    DisjointBoxLayout bl = a_phi.getDBL();
    Box bx = bl.getDomain();
    BoxData<double> phiOut(bx);
    for (int i = 0; i < bl.size();i++)
      {
        a_phi[i].copyTo(phiOut,bl[i]);
      }
    WriteData(phiOut,a_iter,a_dx,a_str.c_str(),a_str.c_str());
  }
} //end namespace Proto
#endif
