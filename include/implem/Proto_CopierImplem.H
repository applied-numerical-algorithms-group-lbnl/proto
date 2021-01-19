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
#if 0
  void FluxRegisterDefine(const DisjointBoxLayout& a_fine,
                          const DisjointBoxLayout& a_coarse,
                          const Point& a_refRatio)
  {
    // Copier to copy flux register increments defined over fine LBD to data 
    // holders defined over the coarse LBD.
    // need assertions here to make sure fine grids are appropriately nested in coarse grids.
    // Loop over Boxes in a_fine to construct fromMotionPlan.
    DataIterator fdit(a_fine);
    int srcProcID = procID();
    for (fdit.begin();!fdit().end();++fdit)
      {
        FineFRIterator ffrit(fdit(),a_refRatio);
        for (ffrit.begin();!ffrit.end();++ffrit)
          {
            auto toIndex = a_coarse.dataIndex(ffrdit());
            auto frBoxSrc = ffrdit.box();
            auto frShift = ffrit.srcShift();
            auto frBoxDest = frBoxSrc.shift(frShift);         
            MotionItem* fromMotionItem = new MotionItem(fdit(),toIndex,frBoxSrc,frBoxDest);
            fromMotionItem->procID = srcProcID;
            fromMotionPlan.push_back(fromMotionItem);
          }
      }
    // Loop over Boxes in a_coarse to construct localMotionPlan, toMotionPlan.
    DataIterator cdit(a_coarse);
    for (cdit.begin();!*cdit != cdit.end();++cdit)
      {
        crseProcID = cdit().procID();
        CoarseFRIterator cfrit(cdit(),a_fine,refRatio);
        for (cfrit.begin();!cfrit.end();++cfrit)
          {
            DataIndex fromIndex = cfrit.fromIndex();
            srcProcID = fromIndex.procID();
            auto frBoxSrc = cfrit.box();
            auto srcShift = cfrit.shift();
            auto frBoxDest = frBoxSrc.shift(srcShift);
            MotionItem* toMotionItem = new MotionItem(fromIndex,cdit(),frBoxSrc,frBoxDest);
            if (crseProcID == srcProcID)
              {
                m_localMotionPlan.push_back(toMotionItem);
              }
            else
              {
                toMotionItem->procID = srcProcID;
                m_toMotionPlan.push_back(toMotionItem);
              }
          }
      }
    sort();
  }
#endif
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
#if 1
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
#endif
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


