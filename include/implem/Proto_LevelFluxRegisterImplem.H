LFRT 
LFR::LevelFluxRegister(const DisjointBoxLayout& a_dblCoarse,
                       const DisjointBoxLayout& a_dblFine,
                       const Point& a_refRatio)
{
  define(a_dblCoarse,a_dblFine,a_refRatio);
}
LFRT 
void LFR::define(const DisjointBoxLayout& a_dblCoarse,
                  const DisjointBoxLayout& a_dblFine,
                  const Point& a_refRatio)
{
  m_refRatio = a_refRatio;
  m_dblFine = a_dblFine;
  m_dblCoarse = a_dblCoarse;
 
  PR_assert(
            m_dblFine.problemDomain().box() 
            == m_dblCoarse.problemDomain().box().refine(m_refRatio));
  // Set up Coarse data holders.
  DataIterator ditC(m_dblCoarse);

  Point finePatchesSize = m_dblCoarse.boxSize()*m_refRatio / m_dblFine.boxSize();
  // Check to see whether Coarse patch is exactly covered by fine patches.
PR_assert(
        finePatchesSize*m_dblFine.boxSize() == m_dblCoarse.boxSize()*m_refRatio);
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
 m_destData.clear();
 m_coarseIndices.clear();
 m_fineRegisters.clear();
  
 int mySizeC = m_dblCoarse.mySize();
 m_coarseRegisters.resize(mySizeC);
 m_destData.resize(mySizeC);
 m_coarseIndices.resize(mySizeC);
 // Set up coarse data holders.  
 for (ditC.begin();*ditC!= ditC.end();++ditC)
   {
     int myindex =m_dblCoarse.myIntIndex(*ditC);
     Point pt = m_dblCoarse[*ditC];
     int k = 0;
     m_coarseRegisters[myindex].clear();
     m_destData[myindex].clear();
     m_coarseIndices[myindex].clear();
     //setup: compute box corresponding to possible fine patches covered by coarse patch.
     Box finePatchPoints = Box(pt,pt).
       refine(m_dblCoarse.boxSize()).refine(m_refRatio).coarsen(m_dblFine.boxSize());
     BoxIterator bit(finePatchPoints);
     // Iterate over the fine patches contained in the coarse box, and see whether patches 
     // adjacent to them in each direction can contribute to a flux register in the coarse 
     // patch. 
     for (bit.begin();bit.done();++bit)
       {
         DataIndex di = m_dblFine.index(*bit);
         for (int d = 0; d < DIM; d++)
           {
             for (int lohi = 0; lohi < 2; lohi++)
               {
                 Point ptAdj = *bit + Point::Basis(d,2*lohi-1);
                 DataIndex diAdj = m_dblFine.find(ptAdj);
                 Box destBox;
                 Box cfBox;
                 if (diAdj != m_dblFine.end() && di == m_dblFine.end())
                   {
                     cfBox = 
                       Box(ptAdj,ptAdj).refine(m_dblFine.boxSize()).coarsen(m_refRatio);
                    
                     if (lohi == 0)
                       {
                         destBox = cfBox.adjCellHi(d,1);
                       }
                     else
                       {
                         destBox = cfBox.adjCellLo(d,1);
                       }
                   }
                 // Intersect with the Problem domain to shift the Box to its proper 
                 // periodic image in the coarse patch.

                 destBox = m_dblCoarse.problemDomain()&destBox; 
                 
                 // Create register and insert it into the coarse register data strucure.
                 
                 int thisKey = key(destBox,di);
                 shared_ptr<BoxData<T,C,MEMTYPE> > temp1(new BoxData<T,C,MEMTYPE>(destBox));
                 shared_ptr<BoxData<T,C,MEMTYPE> > temp2(new BoxData<T,C,MEMTYPE>(destBox));
                 Register<T,C,MEMTYPE> reg(temp1,d,1);
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
 m_hasFineRegister.resize(mySizeF);
 for (int k = 0; k < mySizeF; k++)
   {
     for (int d = 0; d < DIM; d++)
       {
         for (int lohi = 0; lohi < 2; lohi++)
           {
             m_hasFineRegister[k][d][lohi] = -1;
           }
       }
   }
 for (ditF.begin();*ditF!=ditF.end();++ditF)
   {
     Point pt = m_dblFine[*ditF];
     Box bx = Box(pt,pt).refine(m_dblFine.boxSize());
     int myindex = m_dblFine.myIntIndex(*ditF);
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
                 srcBox = bx.adjCellLo(d,1);
                   
               }
             else 
               {
                 srcPt = pt + Point::Basis(d);
                 srcBox = bx.adjCellHi(d,1);
               }
             DataIndex adjIndex = m_dblFine.index(srcPt);
             if (adjIndex == m_dblFine.end())
               {
                 shared_ptr<BoxData<T,C,MEMTYPE> > temp(new BoxData<T,C,MEMTYPE>(srcBox)); 
                 Register<T,C,MEMTYPE> reg(temp,d,lohi);
                 m_fineRegisters[myindex].push_back(reg);
                 m_hasFineRegister[myindex][d][lohi] = k;
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
         Stencil<T>& sten = m_avgFlux[d][lohi] = (0.0)*Shift(Point::Zeros());
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
 // Make sure all the registers are initialized to zero.
 this->reset();
}
LFRT
void LFR::incrementCoarse(
                          const BoxData<T,C,MEMTYPE>& a_flux,
                          const DataIndex& a_di,
                          const T& a_weight,
                          unsigned int a_dir)
{
  int myindex = m_dblCoarse.myIndex(a_di);
  vector<vector<Register<T,C,MEMTYPE > > >& cReg = m_coarseRegisters[myindex];
  Stencil<T> scaleLow  = (-a_weight)*Shift(Point::Zeros() - Point::Basis(a_dir));
  Stencil<T> scaleHigh =   a_weight*Shift(Point::Zeros());
  for (int k = 0; k < cReg.size();k++)
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
                        const BoxData<T,C,MEMTYPE>& a_flux,
                        const DataIndex& a_di,
                        const T& a_weight,
                        unsigned int a_dir)
{
  int myindex = m_dblFine.myIndex(a_di);
  vector<Register<T,C,MEMTYPE> >& myRegs = m_fineRegisters[myindex];
  for (int k = 0; k < myRegs.size(); k++)
    {
      Register<T,C,MEMTYPE>& reg = myRegs[k];
      auto dir = reg.m_dir;
      if (dir == a_dir)
        {
          auto lohi = reg.m_lohi;
          auto db = reg.m_data;
          db += m_avgFlux[a_dir][lohi](a_flux,a_weight);
        }
    }
}
LFRT
void LFR::reflux(LevelBoxData<T,C,MEMTYPE>& a_coarseData,
            const T& a_weight)
{
  // invoke Copier.
  LFRCopierOp<T,C,MEMTYPE> op(*this);
  makeItSo(m_lfrCopier,op);
  // Data has been copied to m_copyDest. We compute the difference between the 
  // coarse and fine registers, and increment a_coarseData.
  Stencil<double> increment = (-a_weight)*Shift(Point::Zeros());
  for (auto dit = a_coarseData.begin();*dit != dit.end();++dit)
    {
      int ind = dit.myIndex();
      int size = m_destData[ind].size();
      vector<shared_ptr<BoxData<T,C> > > dest = m_destData[ind];
      vector<shared_ptr<BoxData<T,C> > > creg = m_coarseRegisters[ind];
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
  for (int myindex = 0; myindex < m_dblCoarse.mySize(); myindex++)
    {
      for (int k = 0; k < m_coarseIndices[myindex].size();k++)
        {
          m_coarseRegisters[myindex][k].m_data->setToZero();
          m_destData[myindex][k]->setToZero();
        }
    }
  for (int myindex =0; myindex < m_dblFine.mySize(); myindex++)
    {
      for (int k = 0; k < m_fineRegisters[myindex].size();k++)
        {
          m_fineRegisters[myindex][k].m_data->setToZero();
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

  Point ptFinePatch = a_bx.low()*m_refRatio/m_dblFine.boxSize();
  Box finePatchPoints = Box(ptFinePatch,ptFinePatch).
    refine(m_dblCoarse.boxSize()).refine(m_refRatio).coarsen(m_dblFine.boxSize());
  Point ptLow = ptFinePatch - Point::Basis(normal);
  Point ptHigh = ptFinePatch + Point::Basis(normal);
  DataIndex diLow = m_dblFine.find(ptLow);
  PR_assert(diLow != m_dblFine.end() ||  m_dblFine.find(ptHigh) != m_dblFine.end());
  int lohi = 1;
  if (diLow != m_dblFine.end()) lohi = 0;
  
  return ind + sz*normal + sz*DIM*lohi;
}
LFRT
BoxData<T,C,MEMTYPE>& LFR::sourceBD(const Box& a_bx,const DataIndex& a_di)
{
  Point ptFine = m_dblFine[a_di];
  Box bxfine = Box(ptFine,ptFine).refine(m_dblFine.boxSize());
  int myindex = m_dblFine.myIndex(a_di);
  int lohi = 2;
  int dir = -1;
  for (int d = 0; d< DIM;d++)
    {
      if (a_bx == bxfine.adjCellLo(d,1))
        {
          lohi = 0;
          dir = d;
        }
      else if  (a_bx == bxfine.adjCellHi(d,1))
        {
          lohi = 1;
          dir = d;
        }
    }
  PR_assert(dir != -1);
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
LFRT
void 
LFR::makeItSo(
              const Copier&       a_copier,
              const LFROP&        a_op) const
{
  int s_verbosity = 0;
  if(s_verbosity > 0) 
    {
      pout() << "makeit so copier = " << endl;
      //a_copier.print();
    }
  {
    PR_TIME("makeItSoBegin");
    makeItSoBegin(a_copier, a_op);
  }
  {
    PR_TIME("makeItSoLocalCopy");
    makeItSoLocalCopy(a_copier, a_op);
  }
  {
    PR_TIME("makeItSoEnd");
    makeItSoEnd(a_op);
  }
}
LFRT
void LFR::makeItSoBegin(
                        const Copier&       a_copier,
                        const LFROP&        a_op) const
{
  // The following five functions are nullOps in uniprocessor mode

#ifdef PR_MPI
  allocateBuffers(a_copier,a_op);  //monkey with buffers, set up 'fromMe' and 'toMe' queues

  writeSendDataFromMeIntoBuffers(a_op);

  // If there is nothing to recv/send, don't go into these functions
  // and allocate memory that will not be freed later.  (ndk)
  // The #ifdef PR_MPI is for the m_buff->m_toMe and m_buff->m_fromMe
  {
    PR_TIME("post_messages");
    this->m_buff->numReceives = m_buff->m_toMe.size();

    if (this->m_buff->numReceives > 0)
      {
        postReceivesToMe(); // all non-blocking
      }  

    this->m_buff->numSends = m_buff->m_fromMe.size();
    if (this->m_buff->numSends > 0)
      {
        postSendsFromMe();  // all non-blocking
      }
  }    
#endif 
}
LFRT
void 
LFR::makeItSoLocalCopy(const Copier&       a_copier,
                       const LFROP&        a_op) const
{
  PR_TIME("local copying");
  CopyIterator it(a_copier, CopyIterator::LOCAL);  
  int items=it.size();
  for (int n=0; n<items; n++)
    {
      const MotionItem& item = it[n];
      //debugging bit to force serial code to run parallel bits
#if 0  
      size_t bufsize_src = a_op.size(item.fromRegion);
      size_t bufsize_dst = a_op.size(item.fromRegion);
      if(bufsize_src != bufsize_dst)
        {
          MayDay<void*>::Error("buffer size mismatch");
        }
      char* charbuffer = new char[bufsize_src];
      a_op.linearOut(charbuffer, item.fromRegion,item.fromIndex);
      a_op.linearIn (charbuffer, item.toRegion,item.toIndex);
        
      delete[] charbuffer;
#else 
      a_op.op(item.fromRegion,item.fromIndex,item.toRegion,item.toIndex);
#endif
        
    }
}
LFRT
void 
LFR::makeItSoEnd(const LFROP& a_op) const
{
#ifdef PR_MPI
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  //cout << "completePendingSends" << " , proc = " << procID() << endl;
  completePendingSends(); // wait for sends from possible previous operation
  //cout << "unpackReceivesToMe"  << " , proc = " << procID() << endl;
  unpackReceivesToMe(a_op); // nullOp in uniprocessor mode
#endif
}
#ifndef PR_MPI
LFRT
void LFR::writeSendDataFromMeIntoBuffers(const LFROP& a_op) const
{
}
LFRT
void LFR::postSendsFromMe() const
{
}
LFRT
void LFR::postReceivesToMe() const
{
}
LFRT
void LFR::unpackReceivesToMe(const LFROP& a_op) const
{
}
#else
//========================================================================
//
// data structures used by makeItSo when we have some
// data that needs to be moved (ie. there are entries
// in the 'FROM' or 'TO' CopyIterators)
//
LFRT
void 
LFR::completePendingSends() const
{
  PR_TIME("completePendingSends");
  //pout() << "numSends " << this->m_buff->numSends << endl;
  if (this->m_buff->numSends > 0)
    {   
      PR_TIME("MPI_Waitall");
      m_buff->m_sendStatus.resize(this->m_buff->numSends);
      //pout() << "this->m_buff->numSends " << this->m_buff->numSends << ", proc = " << procID() << endl;
      //pout() << "&(m_buff->m_sendRequests[0]) " << &(m_buff->m_sendRequests[0]) <<
      // " , " <<  m_buff->m_sendRequests.size() << endl;
        //pout() << "&(m_buff->m_sendStatus[0]) " << &(m_buff->m_sendStatus[0]) <<
      // " , " << m_buff->m_sendStatus.size() << endl; 
      int result = MPI_Waitall(this->m_buff->numSends, &(m_buff->m_sendRequests[0]), &(m_buff->m_sendStatus[0]));
      //pout() << "finished MPI_Waitall" << endl;
      if (result != MPI_SUCCESS)
        {
          cout << "result - failed " << result << endl;
          //hell if I know what to do about failed messaging here
        }
    }
  this->m_buff->numSends = 0;
}
LFRT
void // Getting rid of component dependencies.
LFR::allocateBuffers(
                     const Copier&   a_copier,
                     const LFROP& a_op) const
{
  PR_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  
  if (m_buff->isDefined()) return;
  
  m_buff->m_fromMe.resize(0);
  m_buff->m_toMe.resize(0);
  size_t sendBufferSize = 0;
  size_t recBufferSize  = 0;

  BoxData<T,C> dummy;
  for (CopyIterator it(a_copier, CopyIterator::FROM); it.ok(); ++it)
    {
      const MotionItem& item = it();
      CopierBuffer::bufEntry b;
      b.item = &item;
      b.size = a_op.size(item.fromRegion);
      sendBufferSize+=b.size;
      b.procID = item.procID;
      m_buff->m_fromMe.push_back(b);
    }
  sort(m_buff->m_fromMe.begin(), m_buff->m_fromMe.end());
  for (CopyIterator it(a_copier, CopyIterator::TO); it.ok(); ++it)
    {
      const MotionItem& item = it();
      CopierBuffer::bufEntry b;
      b.item = &item;
      b.size = a_op.size(item.fromRegion);
      recBufferSize+=b.size;
      b.procID = item.procID;
      m_buff->m_toMe.push_back(b);
    }
  sort(m_buff->m_toMe.begin(), m_buff->m_toMe.end());

  // allocate send and receive buffer space.

  if (sendBufferSize > m_buff->m_sendcapacity)
    {
      free((m_buff->m_sendbuffer));
      if (s_verbosity > 0) pout()<<"malloc send buffer "<<sendBufferSize<<std::endl;
      (m_buff->m_sendbuffer) = malloc(sendBufferSize);
      if ((m_buff->m_sendbuffer) == NULL)
        {
          MayDay<void>::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_buff->m_sendcapacity = sendBufferSize;
    }

  if (recBufferSize > m_buff->m_reccapacity)
    {
      free(m_buff->m_recbuffer);
      if (s_verbosity > 0) pout()<<"malloc receive buffer "<<recBufferSize<<std::endl;
      m_buff->m_recbuffer = malloc(recBufferSize);
      if (m_buff->m_recbuffer == NULL)
        {
          MayDay<void>::Error("Out of memory in BoxLayoutData::allocatebuffers");
        }
      m_buff->m_reccapacity = recBufferSize;
    }

  /*
    pout()<<"\n";
    for (int i=0; i<m_buff->m_fromMe.size(); i++)
    pout()<<m_buff->m_fromMe[i].item->region<<"{"<<m_buff->m_fromMe[i].procID<<"}"<<" ";
    pout() <<"::::";
    for (int i=0; i<m_buff->m_toMe.size(); i++)
    pout()<<m_buff->m_toMe[i].item->region<<"{"<<m_buff->m_toMe[i].procID<<"}"<<" ";
    pout() << endl;
  */

  char* nextFree = (char*)(m_buff->m_sendbuffer);
  if (m_buff->m_fromMe.size() > 0)
    {
      for (unsigned int i=0; i<m_buff->m_fromMe.size(); ++i)
        {
          m_buff->m_fromMe[i].bufPtr = nextFree;
          nextFree += m_buff->m_fromMe[i].size;
        }
    }

  nextFree = (char*)m_buff->m_recbuffer;
  if (m_buff->m_toMe.size() > 0)
    {
      for (unsigned int i=0; i<m_buff->m_toMe.size(); ++i)
        {
          m_buff->m_toMe[i].bufPtr = nextFree;
          nextFree += m_buff->m_toMe[i].size;
        }
    }
  // since fromMe and toMe are sorted based on procID, messages can now be grouped
  // together on a per-processor basis.
}
LFRT 
void 
LFR::writeSendDataFromMeIntoBuffers(const LFROP& a_op) const
{
  PR_TIME("write Data to buffers");
  int isize = m_buff->m_fromMe.size();
  for (unsigned int i=0; i< isize; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
      a_op.linearOut(entry.bufPtr,entry.item->fromRegion,entry.item->fromIndex);
    }
}
LFRT
void 
LFR::postSendsFromMe() const
{
  PR_TIME("post_Sends");
  // now we get the magic of message coalescence
  // fromMe has already been sorted in the allocateBuffers() step.
  // PC: Looks like we are just compressing out the buffers by setting the size to zero.
  // This is called after allocateBuffers.
  this->m_buff->numSends = m_buff->m_fromMe.size();
  if (this->m_buff->numSends > 1)
    {
      for (unsigned int i=m_buff->m_fromMe.size()-1; i>0; --i)
        {
          if (m_buff->m_fromMe[i].procID == m_buff->m_fromMe[i-1].procID)
            {
              // Decrementing the number of sends as we coalesce.
              this->m_buff->numSends--;
              m_buff->m_fromMe[i-1].size = m_buff->m_fromMe[i-1].size + m_buff->m_fromMe[i].size;
              m_buff->m_fromMe[i].size = 0;
            }
        }
    }
  m_buff->m_sendRequests.resize(this->m_buff->numSends);
  std::list<MPI_Request> extraRequests;
    
  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->m_buff->numSends; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[next];
      char*  buffer = (char*)entry.bufPtr;
      std::size_t bsize = entry.size;
      int idtag=0;
      while (bsize > PR_MAX_MPI_MESSAGE_SIZE)
        {
          extraRequests.push_back(MPI_Request());
          {
            PR_TIME("MPI_Isend");
            //pout() << "postSendsFromMe -MPI_Isend extra" << endl;
            MPI_Isend(buffer, PR_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                      idtag, Proto_MPI<void>::comm, &(extraRequests.back()));
          }
          maxSize = PR_MAX_MPI_MESSAGE_SIZE;
          bsize -= PR_MAX_MPI_MESSAGE_SIZE;
          buffer+=PR_MAX_MPI_MESSAGE_SIZE;
          idtag++;
        }
      {
        PR_TIME("MPI_Isend");
        //pout() << "postSendsFromMe - MPI_Isend main" << endl;
        MPI_Isend(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Proto_MPI<void>::comm, &(m_buff->m_sendRequests[i]));
      }
      maxSize = std::max<long long>(bsize, maxSize);
      ++next;
      while (next < m_buff->m_fromMe.size() && m_buff->m_fromMe[next].size == 0) ++next;
    }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
    {
      m_buff->m_sendRequests.push_back(*it);
    }
  this->m_buff->numSends = m_buff->m_sendRequests.size();
    
  //PR_MaxMPISendSize = std::max<long long>(PR_MaxMPISendSize, maxSize);
}
LFRT
void 
LFR::postReceivesToMe() const
{
  PR_TIME("post_Receives");
  this->m_buff->numReceives = m_buff->m_toMe.size();
  
  if (this->m_buff->numReceives > 1)
    {
      for (unsigned int i=m_buff->m_toMe.size()-1; i>0; --i)
        {
          if (m_buff->m_toMe[i].procID == m_buff->m_toMe[i-1].procID)
            {
              this->m_buff->numReceives--;
              m_buff->m_toMe[i-1].size += m_buff->m_toMe[i].size;
              m_buff->m_toMe[i].size = 0;
            }
          
        }
    }
  m_buff->m_receiveRequests.resize(this->m_buff->numReceives);
  std::list<MPI_Request> extraRequests;
  unsigned int next=0;
  long long maxSize = 0;
  for (int i=0; i<this->m_buff->numReceives; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_toMe[next];
      char*  buffer = (char*)entry.bufPtr;
      size_t bsize = entry.size;
      int idtag=0;
      while (bsize > PR_MAX_MPI_MESSAGE_SIZE)
        {
          extraRequests.push_back(MPI_Request());
          {
            PR_TIME("MPI_Irecv");
            // pout() << "PostReceivesToMe - MPI_Irecv - extra" << endl;
            MPI_Irecv(buffer,PR_MAX_MPI_MESSAGE_SIZE, MPI_BYTE, entry.procID,
                      idtag, Proto_MPI<void>::comm, &(extraRequests.back()));
          }
          maxSize = PR_MAX_MPI_MESSAGE_SIZE;
          bsize -= PR_MAX_MPI_MESSAGE_SIZE;
          buffer+=PR_MAX_MPI_MESSAGE_SIZE;
          idtag++;
        }
      {
        PR_TIME("MPI_Irecv");
        //pout() << "PostReceivesToMe - MPI_Irecv - main" << endl;
        MPI_Irecv(buffer, bsize, MPI_BYTE, entry.procID,
                  idtag, Proto_MPI<void>::comm, &(m_buff->m_receiveRequests[i]));
      }
      ++next;
      maxSize = std::max<long long>(bsize, maxSize);
      while (next < m_buff->m_toMe.size() && m_buff->m_toMe[next].size == 0) ++next;
    }
  for (std::list<MPI_Request>::iterator it = extraRequests.begin(); it != extraRequests.end(); ++it)
    {
      m_buff->m_receiveRequests.push_back(*it);
    }
  this->m_buff->numReceives = m_buff->m_receiveRequests.size();

  // PR_MaxMPIRecvSize = std::max<long long>(PR_MaxMPIRecvSize, maxSize);
  //pout()<<"maxSize="<<maxSize<<" posted "<<this->m_buff->numReceives<<" receives\n";
}
LFRT
void 
LFR::unpackReceivesToMe(const LFROP& a_op) const
{
  PR_TIME("unpack_messages");
  // pout() << "numReceives " << this->m_buff->numReceives << endl;
  if (this->m_buff->numReceives > 0)
    {
      m_buff->m_receiveStatus.resize(this->m_buff->numReceives);
      int result;
      {
        PR_TIME("MPI_Waitall");
        //pout() << "this->m_buff->numReceives " << this->m_buff->numReceives << endl;
        //pout() << "size (m_buff->m_receiveRequests[0]) " << &(m_buff->m_receiveRequests[0]) <<
        //  " , " <<  m_buff->m_receiveRequests.size() << endl;
        //pout() << "&(m_buff->m_receiveStatus[0]) " << &(m_buff->m_receiveStatus[0]) <<
        //  " , " << m_buff->m_receiveStatus.size() << endl; 
        result = MPI_Waitall(this->m_buff->numReceives, &(m_buff->m_receiveRequests[0]),
                             &(m_buff->m_receiveStatus[0]));
        //pout() << "unpackReceives - post Waitall " << this->m_buff->numReceives << endl;
      }
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
          //maybe a mayday::warning?
        }

      int isize = m_buff->m_toMe.size();
      for (unsigned int i=0; i< isize; ++i)
        {
          const CopierBuffer::bufEntry& entry = m_buff->m_toMe[i];
          a_op.linearIn(entry.bufPtr,entry.item->toRegion,entry.item->toIndex);
        }
    }
  this->m_buff->numReceives = 0;
}
#endif
