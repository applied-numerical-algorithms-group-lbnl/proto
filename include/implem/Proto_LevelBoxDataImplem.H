#define LBDT template \
  <typename T, unsigned int C, MemType MEMTYPE, unsigned int CENTERING>
LBDT
LBD::LevelBoxData()
{
  m_isDefined = false;
};
///
LBDT
LBD::
LevelBoxData(const DisjointBoxLayout& a_dbl, 
             const Point& a_grow)
{
  define(a_dbl, a_grow);
};
LBDT
void LBD::define(
                 const DisjointBoxLayout& a_dbl, 
                 const Point& a_ghost)
{
  m_isDefined = true;
  m_ghost = a_ghost;
  m_dbl = a_dbl;
  DataIterator dit(m_dbl);
  for (dit.begin();*dit != dit.end();++dit)
    {
      shared_ptr<BoxData<T,C> > temp(new BoxData<T,C>(dit.box().grow(m_ghost)));
      m_data.push_back(temp);
    }
  
  if (m_ghost != Point::Zeros())
    {
      m_exchangeCopier.exchangeDefine(m_dbl,m_ghost);
      // pout() << m_exchangeCopier << endl;
    }
};
/// PC::implement
LBDT
unsigned int 
LBD::size() const
{
  return m_data.size();
}

///PC::implement
LBDT
void   
LBD::setToZero()
{
  for(unsigned int ibox = 0; ibox < m_data.size(); ibox++)
    {
      m_data[ibox]->setVal(0);
    }
};
/// fills ghost cells with valid data
// PC:: need to figure out LDOperator.
LBDT
void 
LBD::exchange()
{
  PR_TIME("exchange");
#ifdef PR_MPI
  {
      PR_TIME("MPI_Barrier exchange");
      MPI_Barrier(Proto_MPI<void>::comm);
  }
#endif
  this->makeItSo(*this, *this, m_exchangeCopier);
}
/** call exchange if you want ghost cells filled */
LBDT
void 
LBD::
copyTo(LBD& a_dest) const
{
  Copier copytocopier;
  copytocopier.LBDCopyToDefine(m_dbl,a_dest.m_dbl);
  //pout() << copytocopier << endl;
  this->makeItSo(*this, a_dest, copytocopier);
}
//// Interaction with iterators.
LBDT
DataIterator LBD::begin() const
{
  DataIterator retval(m_dbl);
  retval.begin();
  return retval;
} 
LBDT
BoxData<T,C>& 
LBD::operator[](const DataIndex& a_di) 
{
  return (*(m_data[m_dbl.myIntIndex(a_di)]));
}
LBDT
const BoxData<T,C>& 
LBD::operator[](const DataIndex& a_di) const
{
  return (*(m_data[m_dbl.myIntIndex(a_di)]));
}
LBDT
void 
LBD::makeItSo(const LBD&          a_src,
              LBD&                a_dest,
              const Copier&       a_copier,
              const LDOP&         a_op) const
{
  if(s_verbosity > 0) 
    {
      pout() << "makeit so copier = " << endl;
      //a_copier.print();
    }
  {
    PR_TIME("makeItSoBegin");
    makeItSoBegin(a_src, a_dest, a_copier, a_op);
  }
  {
    PR_TIME("makeItSoLocalCopy");
    makeItSoLocalCopy(a_src, a_dest, a_copier, a_op);
  }
  {
    PR_TIME("makeItSoEnd");
    a_dest.makeItSoEnd(a_op);
  }
}
LBDT
void LBD::makeItSoBegin(const LBD&          a_src,
                        LBD&                a_dest,
                        const Copier&       a_copier,
                        const LDOP&         a_op) const
{
  // The following five functions are nullOps in uniprocessor mode

#ifdef PR_MPI
  allocateBuffers(a_src,a_dest,a_copier,a_op);  //monkey with buffers, set up 'fromMe' and 'toMe' queues

  writeSendDataFromMeIntoBuffers(a_src,a_op);

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
LBDT
void 
LBD::makeItSoLocalCopy(const LBD&          a_src,
                       LBD&                a_dest,
                       const Copier&       a_copier,
                       const LDOP&         a_op) const
{
  PR_TIME("local copying");
  CopyIterator it(a_copier, CopyIterator::LOCAL);  
  int items=it.size();
  for (int n=0; n<items; n++)
    {
      const MotionItem& item = it[n];
      //debugging bit to force serial code to run parallel bits
      // PC:: check with BVS
#if 0
      const T & srcFAB = a_src[item.fromIndex];
      T       & dstFAB = a_dest[item.toIndex];
      size_t bufsize_src = a_op.size(srcFAB, item.fromRegion);
      size_t bufsize_dst = a_op.size(srcFAB, item.fromRegion);
      if(bufsize_src != bufsize_dst)
        {
          MayDay::Error("buffer size mismatch");
        }
      char* charbuffer = new char[bufsize_src];
      a_op.linearOut(srcFAB, charbuffer, item.fromRegion);
      a_op.linearIn (dstFAB, charbuffer, item.toRegion);
        
      delete[] charbuffer;
#else
      a_op.op(a_dest[item.toIndex], 
              item.fromRegion,
              item.toRegion,
              a_src[item.fromIndex]);
#endif
        
    }
}
LBDT
void 
LBD::makeItSoEnd(const LDOP& a_op)
{
  // Uncomment and Move this out of unpackReceivesToMe()  (ndk)
  //cout << "completePendingSends" << " , proc = " << procID() << endl;
  completePendingSends(); // wait for sends from possible previous operation
  //cout << "unpackReceivesToMe"  << " , proc = " << procID() << endl;
  unpackReceivesToMe(a_op); // nullOp in uniprocessor mode
}
#ifndef PR_MPI
// uniprocessor version of all these nullop functions.
LBDT
void LBD::completePendingSends() const
{
  PR_TIME("completePendingSends");
  if (this->m_buff->numSends > 0)
    {
      PR_TIME("MPI_Waitall");
      m_buff->m_sendStatus.resize(this->m_buff->numSends);
      int result = MPI_Waitall(this->m_buff->numSends, &(m_buff->m_sendRequests[0]), &(m_buff->m_sendStatus[0]));
      if (result != MPI_SUCCESS)
        {
          //hell if I know what to do about failed messaging here
        }
    }
  this->m_buff->numSends = 0;
}
// Need to get rid of components.
#if 0
LBDT
void LBD::allocateBuffers(const LBD& a_src, 
                          const LBD& a_dest,
                          const Copier& a_copier,
                          const LDOP& a_op) const
{
  PR_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  a_dest.m_buff = m_buff;
  
  if (m_buff->isDefined()) return;

  if(s_verbosity > 0)
    {
      pout() << " allocate buffers "  << endl;
    }

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
      b.size = a_op.size(a_src[item.fromIndex], item.fromRegion);
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
      b.size = a_op.size(a_dest[item.toIndex], item.fromRegion);
      recBufferSize+=b.size;
      b.procID = item.procID;
      m_buff->m_toMe.push_back(b);
    }
  sort(m_buff->m_toMe.begin(), m_buff->m_toMe.end());

  // allocate send and receveive buffer space.

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
#endif
LBDT
void LBD::writeSendDataFromMeIntoBuffers(const LBD& a_src,foo
                                         const LDOP& a_op) const
{
}

LBDT
void LBD::postSendsFromMe() constfoo
{
}

LBDT
void LBD::postReceivesToMe() constfoo
{
}
LBDT
void LBD::unpackReceivesToMe(const LDOP& a_op) foo
{
}
#else
//========================================================================
//
// data structures used by makeItSo when we have some
// data that needs to be moved (ie. there are entries
// in the 'FROM' or 'TO' CopyIterators)
//
LBDT
void 
LBD::completePendingSends() const
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
LBDT
void // Getting rid of component dependencies.
LBD::allocateBuffers(const LBD& a_src,
                     const LBD& a_dest,
                     const Copier&   a_copier,
                     const LDOP& a_op) const
{
  PR_TIME("MPI_allocateBuffers");
  m_buff = &(((Copier&)a_copier).m_buffers);
  a_dest.m_buff = m_buff;
  
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
      b.size = a_op.size(item.fromRegion); // ?? Need to look at a_op API.
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
LBDT 
void 
LBD::writeSendDataFromMeIntoBuffers(const LBD& a_src, 
                                    const LDOP& a_op) const
{
  PR_TIME("write Data to buffers");
  int isize = m_buff->m_fromMe.size();
  for (unsigned int i=0; i< isize; ++i)
    {
      const CopierBuffer::bufEntry& entry = m_buff->m_fromMe[i];
      a_op.linearOut(a_src[entry.item->fromIndex], entry.bufPtr,
                     entry.item->fromRegion);
    }
}
LBDT
void 
LBD::postSendsFromMe() const
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
LBDT
void 
LBD::postReceivesToMe() const
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
LBDT
void 
LBD::unpackReceivesToMe(const LDOP& a_op) 
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
          a_op.linearIn(this->operator[](entry.item->toIndex), entry.bufPtr,
                        entry.item->toRegion);
        }
    }
  this->m_buff->numReceives = 0;
}
#endif
#if 0
void WriteData(LevelBoxData<double>& a_phi,int a_iter,double a_dx,std::string a_str)
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
#endif
