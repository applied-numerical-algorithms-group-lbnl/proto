#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "Proto_Copier.H"
#include <chrono>
#include <vector>

using std::ostream;

namespace Proto
{
  CopierBuffer::~CopierBuffer()
  {
    clear();
  }

  void CopierBuffer::clear()
  {
    if (m_sendbuffer != NULL) free(m_sendbuffer);
    if (m_recbuffer  != NULL) free(m_recbuffer);
    m_sendbuffer = NULL;
    m_recbuffer  = NULL;
    m_sendcapacity = 0;
    m_reccapacity = 0;
    m_isDefined = false;
  }
  bool Copier::operator==(const Copier& rhs) const
  {
    if(m_localMotionPlan.size() != rhs.m_localMotionPlan.size()) return false;
    if(m_fromMotionPlan.size() != rhs.m_fromMotionPlan.size()) return false;
    if(m_toMotionPlan.size() != rhs.m_toMotionPlan.size()) return false;
    for(int i=0; i<m_localMotionPlan.size(); ++i)
      {
        if(!(*m_localMotionPlan[i] == *rhs.m_localMotionPlan[i]))
          {
            return false;
          }
      }
    return true;
  }

  void Copier::LBDCopyToDefine(const DisjointBoxLayout& a_src,
                               const DisjointBoxLayout& a_dest)
  {
    // Loop over Boxes in the destination to construct toMotionPlan, localMotionPlan. 
    DataIterator ddit(a_dest);
    int destProcID = procID();
    for (ddit.begin();*ddit != ddit.end();++ddit)
      {
        Box bx = ddit.box();
        NeighborIterator nit(a_src,bx);
        for (nit.begin();*nit != nit.end();++nit)
          {
            int nbrProcID = nit.procid();
            Box destBx = nit.destBox();
              if (!destBx.empty())
                { 
                  MotionItem* toMotionItem = new MotionItem(*nit, *ddit, destBx, destBx);
                  if (nbrProcID == procID())
                    { // local move
                      m_localMotionPlan.push_back(toMotionItem);
                    }
                  else
                    {
                      toMotionItem->procID = nbrProcID;
                      m_toMotionPlan.push_back(toMotionItem);
                    }
                }
          }
      }
    // Iterate over source DBL to obtain entries in m_fromMotionPlan.           
    DataIterator sdit(a_src);
    //  unsigned int srcProcID = procID();
    for (sdit.begin();*sdit != sdit.end();++sdit)
      {
        Box bx = sdit.box();
        NeighborIterator nit(a_dest,bx);
        for (nit.begin();*nit != nit.end();++nit)
          {
            int nbrProcID = nit.procid();
            Box srcBx = nit.srcBox();
            if (nbrProcID != destProcID)
              {
                MotionItem* fromMotionItem = new MotionItem(*sdit,*nit,srcBx,srcBx);
                fromMotionItem->procID = nbrProcID;
                m_fromMotionPlan.push_back(fromMotionItem);
              }
          }
      }     
    sort();
  }
  void Copier::exchangeDefine(const DisjointBoxLayout& a_src, 
                              const Point&  a_ghost)
  {
    // Loop over all Boxes in the DBL.
    DataIterator dit(a_src);
    int myProcID = procID();
    for (dit.begin();*dit != dit.end();++dit)
      {
        Box bx0 = dit.box();
        Box bx = a_src.problemDomain() & dit.box().grow(a_ghost);
        NeighborIterator nit(a_src,bx);
        for (nit.begin();*nit!=nit.end();++nit)
          {
            int nbrProcID = nit.procid();
            Point nbrShift = nit.shift();
            // For each pair of interacting boxes, we construct to, from and local. 
            // Compute the boxes corresponding to filling the ghost region from *dit to *nit.
            // The source box is shifted, the dest box is not. 
            {  
              Box destBx = nit.destBox();
              if (!destBx.empty() && !bx0.contains(destBx))
                {
                  Box srcBx = nit.srcBox(); 
                  MotionItem* toMotionItem = new MotionItem(*nit, *dit, srcBx, destBx);
                  if (nbrProcID == myProcID)
                    { // local move
                      m_localMotionPlan.push_back(toMotionItem);
                    }
                  else
                    {
                      toMotionItem->procID = nbrProcID;
                      m_toMotionPlan.push_back(toMotionItem);
                    }
                }
            }
            // Construct Motion item for filling the ghost region for *nit from *dit.       
            if (nbrProcID != myProcID)
              {
                Box nbrbx = Box(nit.point(),nit.point()).refine(a_src.boxSize()).grow(a_ghost);
                Box nbrbxShift = nbrbx.shift(-nbrShift);       
                Box srcFromBx = nbrbxShift&bx0;
                Box destFromBx = srcFromBx.shift(nbrShift);
                MotionItem* fromMotionItem = new MotionItem(*dit, *nit, srcFromBx, destFromBx);
                fromMotionItem->procID = nbrProcID;
                m_fromMotionPlan.push_back(fromMotionItem);
              }
          } 
      }   
    sort();
  }
 template<typename T,unsigned int C, MemType MEMTYPE>
 void 
 Copier::LFRCopierDefine(const DisjointBoxLayout& a_dblFine,
                         const DisjointBoxLayout& a_dblCoarse,
                         const vector<vector<Register<T,C,MEMTYPE > > >& a_coarseRegisters,
                         const vector<array<array<int,2>, DIM> >& a_hasFineRegister,
                         const Point& a_refRatio)
  {
    // Copier to copy flux register increments defined over fine LBD to data 
    // holders defined over the coarse LBD.

    // need assertions here to make sure fine grids are appropriately nested in coarse grids.

    // Loop over coarse Boxes I own to construct toMotionPlan (ones I am recieving)
    // and localMotionPlan.
    DataIterator ddit(a_dblCoarse);
    for (ddit.begin();*ddit != ddit.end();++ddit)
      {
        int myindex = a_dblCoarse.myIntIndex(*ddit);
        //const vector<vector<Register<T,C,MEMTYPE > > >& cReg = a_coarseRegisters[myindex];
        auto& cReg = a_coarseRegisters[myindex];
        for (int k = 0; k < cReg.size(); k++)
          {
            const Register<T,C,MEMTYPE >& reg = cReg[k];
            // Find fine iterator corresponding to this flux register.
            Box destBox = reg.m_data->box();
            int lohi = reg.m_lohi;
            int dir = reg.m_dir;
            Point ptFine = destBox.low()*a_refRatio/a_dblFine.boxSize() + Point::Basis(dir,1-2*lohi);
            DataIndex diFine = a_dblFine.find(ptFine);
            PR_assert(diFine != a_dblFine.end());
            // Compute fine register box (nontrivial b/c of possibility of periodic BCs. 
            Box coarsenedFineBox = a_dblFine.box(diFine).coarsen(a_refRatio);
            Box srcBox;
            if (lohi == 0)
              {
                srcBox = coarsenedFineBox.adjCellLo(dir,1);
              }
            else
              {
                srcBox =  coarsenedFineBox.adjCellHi(dir,1);
              }                       
            MotionItem* toMotionItem = new MotionItem(diFine, *ddit, srcBox, destBox);
            if (a_dblFine.procid(diFine) == procID())
              { 
                m_localMotionPlan.push_back(toMotionItem);
              }
            else
              {
                toMotionItem->procID =a_dblFine.procid(diFine);
                m_toMotionPlan.push_back(toMotionItem);
              }
          }
      }
    // Iterate over fine grids that I own to obtain entries in m_fromMotionPlan.           
    DataIterator sdit(a_dblFine);
    int destProcID = procID();
    for (sdit.begin();*sdit != sdit.end();++sdit)
      {
        int myindex = a_dblFine.myIntIndex(*sdit);
        Box coarsenedFineBox = a_dblFine.box(*sdit).coarsen(a_refRatio);
        for (int dir = 0; dir < DIM; dir++)
          {
            for (int lohi = 0;lohi < 2; lohi++)
              {
                if (a_hasFineRegister[myindex][dir][lohi] < 2)
                  {
                    Box srcBox;
                    // compute box for possible register.
                    if (lohi == 0)
                      {
                        srcBox = coarsenedFineBox.adjCellLo(dir,1);
                      }
                    else
                      {
                        srcBox = coarsenedFineBox.adjCellHi(dir,1);
                      }
                    auto ddit = a_dblCoarse.find(srcBox.coarsen(a_dblCoarse.boxSize()).low());                   
                    int coarseProcID = a_dblCoarse.procid(ddit);
                    Box destBox =  
                      srcBox.shift(a_dblCoarse.problemDomain().shifted(srcBox.low()));
                    if (coarseProcID != procID())
                      {
                        MotionItem* fromMotionItem = 
                          new MotionItem(*sdit,ddit,srcBox,destBox);
                        fromMotionItem->procID = coarseProcID;
                        m_fromMotionPlan.push_back(fromMotionItem);
                      }
                  }
              }
          }
      }     
    sort();
  }
  int Copier::print() const
  {
    pout()  << *this;
    return 0;
  }

  int Copier::numLocalCellsToCopy() const
  {
    int sum = 0;
    for (unsigned int i = 0; i < m_localMotionPlan.size(); ++i)
      {
        sum += m_localMotionPlan[i]->fromRegion.size();
      }
    return sum;
  }

  int Copier::numFromCellsToCopy() const
  {
    int sum = 0;
    for (unsigned int i = 0; i < m_fromMotionPlan.size(); ++i)
      {
        sum += m_fromMotionPlan[i]->fromRegion.size();
      }
    return sum;
  }

  int Copier::numToCellsToCopy() const
  {
    int sum = 0;
    for (unsigned int i = 0; i < m_toMotionPlan.size(); ++i)
      {
        sum += m_toMotionPlan[i]->fromRegion.size();
      }
    return sum;
  }

  ostream& operator<< (ostream& os, const Copier& copier)
  {
    os << "local(" << procID() << "): ";
    for (CopyIterator it(copier, CopyIterator::LOCAL); it.ok(); ++it)
      {
        os << " from " << it().fromIndex << it().fromRegion << " to " << it().toIndex << it().toRegion << '\n'
           << "          ";
      }
    os << "\nfrom(" << procID() << "): ";
    for (CopyIterator it(copier, CopyIterator::FROM); it.ok(); ++it)
      {
        os << " from " << it().fromIndex << it().fromRegion << " to " << it().toIndex << it().toRegion << "[" << it().procID << "]" << '\n'
           << "         ";
      }
    os << "\nto(" << procID() << "): ";
    for (CopyIterator it(copier, CopyIterator::TO); it.ok(); ++it)
      {
        os << " from " << it().fromIndex << it().fromRegion << " to " << it().toIndex << it().toRegion << "[" <<it().procID << "]" << '\n'
           << "       ";
      }
    os << "\n";
    return os;
  }
  Copier::~Copier()
  {
    PR_TIME("~Copier");
    clear();
  }

  void Copier::setBufferAllocated(bool arg) const
  {
    buffersAllocated  = arg;
  }

  bool Copier::bufferAllocated() const
  {
    return buffersAllocated;
  }
  class MotionItemSorter
  {
  public:
    inline bool operator()(MotionItem* const & lhs, MotionItem* const& rhs) const;
  };
  
  inline bool  MotionItemSorter::operator()(MotionItem* const & lhs, MotionItem* const& rhs) const
  {
    if (lhs->procID == rhs->procID)
      {
        const Box& left = lhs->toRegion;
        const Box& right = rhs->toRegion;
        if (left.low() == right.low())
          {
            return left.high() < (right.high()); //PC:: LexLT??
          }
        else
          {
            return left < right;
          }
      }
    return lhs->procID < rhs->procID;
  }
  void Copier::sort()
  {
    std::vector<MotionItem*>& vfrom  = m_fromMotionPlan;
    std::sort(vfrom.begin(), vfrom.end(), MotionItemSorter());
    std::vector<MotionItem*>& vto = m_toMotionPlan;
    std::sort(vto.begin(), vto.end(), MotionItemSorter());
  }
void Copier::clear()
{
  for (unsigned int i = 0; i < m_localMotionPlan.size(); ++i)
  {
    delete m_localMotionPlan[i];
  }
  for (unsigned int i = 0; i < m_fromMotionPlan.size(); ++i)
  {
    delete m_fromMotionPlan[i];
  }
  for (unsigned int i = 0; i < m_toMotionPlan.size(); ++i)
  {
    delete m_toMotionPlan[i];
  }

  m_localMotionPlan.resize(0);
  m_fromMotionPlan.resize(0);
  m_toMotionPlan.resize(0);
  m_isDefined = false;
  m_buffers.clear();
}
}
