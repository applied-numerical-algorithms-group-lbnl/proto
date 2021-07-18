LFRT 
LFR::LevelFluxRegister(const DisjointBoxLayout& a_dblFine,
                       const DisjointBoxLayout& a_dblCoarse,
                       const Point& a_refRatio)
{
  define(a_dblFine,a_dblCoarse,a_refRatio);
}
LFRT 
void LFR::define(const DisjointBoxLayout& a_dblFine,
                  const DisjointBoxLayout& a_dblCoarse,
                  const Point& a_refRatio)
{
  m_refRatio = a_refRatio;
  m_dblFine = a_dblFine;
  m_dblCoarse = a_dblCoarse;
 
  
  // Set up Coarse data holders.
  DataIterator ditC(m_dblCoarse);

  m_finePatchesBoxSize = m_dblCoarse.boxsize()*m_refRatio / m_dblFine.boxsize();
  // Check to see whether Coarse patch is exactly covered by fine patches.
PR_assert(
        m_finePatchesSize*m_dblFine.boxSize() == m_dblCoarse.boxSize()*m_refRatio);
// Check to see whether coarsened fine patch is more than one coarse cell wide in all directions that are being refined. 

// Fine registers are defined relative to the fine patch coordinates, and coarse 
// registers are defined relative to the coarse patch coordinates. What happens in 
// the presence of periodic boundaries is mediated in LFROperator.

 for (int dir = 0; dir < DIM; dir++)
   {
     if (m_refRatio[dir] > 1)
       {
         PR_assert(m_dblFine.boxSize()[dir]/m_refRatio[dir] > 1);
       }
   }
 m_coarseRegisters.clear();
 m_copyDest.clear();
 m_coarseIndices.clear();
 m_fineRegisters.clear();
 m_fineIndices.clear();
  
 int mySizeC = m_dblCoarse.mySize();
 m_coarseRegisters.resize(mySizeC);
 m_copyDest.resize(mySizeC);
 m_coarseIndices.resize(mySizeC);
 // Set up coarse data holders.  
 for (ditC.begin();*ditC!= ditC.end();++ditC)
   {
     int myindex = ditC.myIndex();
     Point pt = m_dblCoarse[*ditC];
     int k = 0;
     m_coarseRegisters[myindex].clear();
     m_copyDest[myIndex].clear();
     m_coarseIndices[myIndex].clear();
     //setup: compute box corresponding to possible fine patches covered by coarse patch.
     Box finePatchPoints = Box(pt,pt).
       refine(m_dblCoarse.boxSize()).refine(m_refRatio).coarsen(m_dblFine.boxSize());
     BoxIterator bit(finePatchPoints);
     // Iterate over the fine patches contained in the coarse box, and see whether patches 
     // adjacent to them in each direction can contribute to a flux register in the coarse 
     // patch. 
     for (bit.begin();bit.done();++bit)
       {
         DataIndex di = m_dblFine(*bit);
         for (int d = 0; d < DIM; d++)
           {
             for (int lohi = 0; lohi < 2;lohi)
               {
                 Point ptAdj = *bit + Point::Basis(d,2*lohi-1);
                 DataIndex diAdj = m_dblFine.find(ptLow);
                 if (diAdj != m_dblFine.end() && di == m_dblFine.end())
                   {
                     Box cfBox = 
                       Box(ptAdj,ptAdj).refine(m_dblFine.boxSize()).coarsen(m_refRatio);
                     Box destBox;
                     if (lohi = 0)
                       {
                         destBox = cfBox.adjCellHi(d);
                       }
                     else
                       {
                         destBox = cfBox.adjCellLo(d);
                       }
                   }
                 // Intersect with the Problem domain to shift the Box to its proper 
                 // periodic image in the coarse patch.

                 destBox = m_coarseDBL.problemDomain()&destBox; 
                 
                 // Create register and insert it into the coarse register data strucure.
                 
                 int thisKey = key(cfBox,di);
                 shared_ptr<BoxData<T,C> > temp1(new BoxData<T,C>(destBox));
                 shared_ptr<BoxData<T,C> > temp2(new BoxData<T,C>(destBox));
                 Register reg(temp1,d,1);
                 m_coarseRegisters[myindex].push_back(reg);
                 m_destData[myindex].push_back(temp2);
                 m_coarseIndices[myindex][thisKey] = k;
                 k++;
               }
           }
       }
   }  
 
  // Set up fine data holders.
 DataIterator ditF(m_dblFine);
 int mySizeF = m_dblFine.mySize();
 for (ditF.begin();*ditF!=ditF.end();++ditF)
   {
     Point pt = m_dblFine[*ditF];
     Box bx = Box(pt,pt).refine(m_dblFine.boxSize());
     int myIndex = m_dblFine.myIndex(*ditF);
     int k = 0;
     Box srcBox;
     for (int d = 0; d < DIM; d++)
       {
         for (int lohi=0;lohi < 2; lohi++)
           {
             Point srcPt;
             if (lohi == 0) 
               {
                 srcPt = pt + Point::Basis(d,-1);
                 srcBox = bx.adjCellLo(d);
                   
               }
             else 
               {
                 srcPt = pt + Point::Basis(d);
                 srcBox = bx.adjCellHi(d);
               }
             DataIndex adjIndex = m_dblFine(srcPt);
             if (adjIndex == m_dblFine.end())
               {
                 shared_ptr<BoxData<T,C> > temp(new BoxData<T,C>(srcBox)); 
                 Register reg(temp,d,lohi);
                 m_fineRegisters[myIndex].push_back(reg);
                 m_hasRegisters[myIndex][d][lohi] = k;
                 k++;
               }
           }
       }
   }

// Now that we have created the Register data holders, we can define the copier.
 m_lfrCopier.LFRCopierDefine(m_dblFine,m_dblCoarse,
                             m_coarseRegisters,m_hasFineRegister,m_refRatio);

 // Create stencils for averaging down fluxes.
 for (int d = 0; d < DIM; d++)
   {
     for (int lohi = 0; lohi < 1; lohi++)
       {
         Stencil& sten = m_avgFlux[d][lohi] = (0.0)*Shift(Point::Zeros());
         Box bx(Point::Zeros(),(Point::Ones() - Point::Basis(d))*m_refRatio);
         BoxIterator bit(bx);
         T coef = 1.0/bx.size();
         for (auto iter = bx.begin(); iter !=bx.end(); ++iter)
           {
             sten += coef*Shift(*iter);
           }
         sten.srcRatio() = m_refRatio; 
         if (lohi == 0)
           {
             sten.destShift() = Point::Basis(d,-1);
           }
       }
   }
}
LFRT
void LFR::incrementCoarse(const BoxData<T,C>& a_flux,
                          const DataIndex& a_di,
                          const T& a_weight,
                          int a_dir)
{
  int myIndex = a_di.myIndex();
  vector<vector<Register> >& cReg = m_coarseRegisters[myIndex];
  Stencil scaleLow  = (-a_weight)*Shift(Point::Zeros() - Point::Basis(a_dir));
  Stencil scaleHigh =   a_weight*Shift(Point::Zeros());
  for (int k = 0; k < cRegs.size();k++)
    {   
      if (cReg[k].m_dir == a_dir)
        {
          BoxData<T,C>& reg = cReg[k].m_data; 
          int side = cReg[k].m_lohi;
          if (side == 0)
            {
              reg += scaleLow(a_flux);
            }
          else
            {
              reg += scaleHigh(a_flux);
            }
        }
    }
}
LFRT
void LFR::incrementFine(
                        const BoxData<T,C>& a_flux,
                        const DataIndex& a_di,
                        const T& a_weight,
                        int a_dir)
{
  int myIndex = a_di.myIndex();
  vector<BoxData<T,C> >& myRegs = m_fineRegisters[myIndex];
  for (int dir = 0; k < myRegs.size(); k++)
    {     
      auto dir = register.m_dir;
      if (m_dir == a_dir)
        {
          auto register = myRegs[k];
          auto lohi = register.m_lohi;
          auto db = register.m_data;
          db += m_avgFlux[d][lohi](a_flux);
        }
    }
}
LFRT
void LFR::reflux(LevelBoxData<T,C>& a_coarseData,
            const T& a_weight)
{
  // invoke Copier.
  LFRCopierOp op(*this);
  makeItSo<LFRCopier>(m_lfrCopier,op);
  // Data has been copied to m_copyDest. We compute the difference between the 
  // coarse and fine registers, and increment a_coarseData.
  Stencil increment = (-a_weight)*Shift(Point::Zeros());
  for (dit = a_coarseData.begin();*dit != dit.end();++dit)
    {
      int ind = dit.myIndex();
      int size = m_coarseDest[ind].size();
      vector<shared_ptr<BoxData<T,C> > dest = m_copyDest[ind];
      vector<shared_ptr<BoxData<T,C> > creg = m_coarseRegisters[ind];
      for (int k = 0; k < size; k++)
        {
          *(dest[k]) += increment(*(creg[k]));  
        }
      for (int k = 0; k < size; k++)
        {
          a_coarseData[*dit] += *(creg[k]);
        }
    }          
}
LFRT
void LFR::reset()
{
  // set values in the registers to zero.
  for (dit = m_dblCoarse.begin();*dit != dit.end(); ++dit)
    {
      int myIndex = dit.myIndex();
      for (int k = 0; k < m_coarseIndices[myIndex].size();k++)
        {
          m_coarseRegisters[myIndex][k].m_data->setToZero();
          m_copyDest[myIndex][k]->setToZero();
        }
    }
  for (dit = m_dblFine.begin();*dit != dit.end(); ++dit)
    {
      int myIndex = dit.myIndex();
      for (int k = 0; k < m_fineIndices[myIndex].size();k++)
        {
          m_fineRegisters[myIndex][k].m_data->setToZero();
        }
    }
}
LFRT
int LFR::key(const Box& a_bx,const DataIndex& a_di)
{
  // Find which direction is normal to the register. We use the fact that the register
  // is of size 1 only in the normal direction.

  int normal = -1;
  for (int dir = 0; dir < DIM; dir++)
    {
      if (a_bx.size(dir) == 1) normal = dir;
    }
  PR_assert(normal !=-1);
  Point coarsePt = m_dblCoarse[a_di];
  Box bxCoarse(Point::Zeros(),m_dblCoarse.boxSize());
  Point shiftedPt = a_bx.low() - coarsePt;
  unsigned int ind = bxCoarse.index(shiftedPt);
  unsigned int sz = bxCoarse.size();

  // Find whether this is a low or high side register. We use the fact that, if 
  // fine grids are at least two coarse cells wide, there can be only one fine patch 
  // per direction adjacent to the register Box.

  Point ptFinePatch = a_bx.low()*m_refine/m_dblFine.boxSize();
  Box finePatchPoints = Box(pt,pt).
    refine(m_dblCoarse.boxSize()).refine(m_refRatio).coarsen(m_dblFine.boxSize());
  Point ptLow = ptFinePatch - Point::Basis(normal);
  Point ptHigh = ptFinePatch + Point::Basis(normal);
  DataIndex diLow = m_dblFine.find(ptLow);
  PR_assert(diLow != m_dblFine.end() ||  m_dblFine.find(ptHigh) != dblFine.end());
  lohi = 1;
  if (diLow != m_dblFine.end()) lohi = 0;
  
  return ind + sz*normal + sz*DIM*lohi;
}
LFRT
BoxData<T,C,MEMTYPE>& LFR::sourceBD(const Box& a_bx,const DataIndex& a_di)
{
  Point ptFine = m_dblFine[a_di];
  Box bxfine = Box(ptFine,ptFine).refine(m_dblFine.boxSize());
  int myindex = m_dblFine.myIndex(a_di);
  lohi = -1;
  dir = -1;
  for (int d = 0; d< DIM;d++)
    {
      if (a_bx == bxFine.adjCellLo(d))
        {
          lohi = 0;
          dir = d;
        }
      else if  (a_bx == bxFine.adjCellHi(d))
        {
          lohi = 1;
          dir = d;
        }
    }
  PR_ASSERT(dir != -1);
  int k = m_hasFineRegister[dir][lohi];
  return *(m_fineRegisters[myindex][k].m_data);
}
LFRT
BoxData<T,C,MEMTYPE>& LFR::destBD(const Box& a_bx,const DataIndex& a_di)
{
  int keyCoarse = key(a_bx,a_di);
  int myindex = m_dblCoarse.myIndex(a_di);
  return *(m_coarseRegisters[myindex][keyCoarse].m_data);
}
