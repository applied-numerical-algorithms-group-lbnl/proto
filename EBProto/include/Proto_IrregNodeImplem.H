namespace Proto
{
  ////
  int IrregNode::
  index(int a_idir, Side::LoHiSide a_sd)
  {
    CH_assert(a_idir >= 0 && a_idir < DIM);
    int retval;
    if (a_sd == Side::Lo)
    {
      retval = a_idir;
    }
    else
    {
      retval = a_idir + DIM;
    }
    return retval;
  }

  ///
  void 
  IrregNode::
  setMomentsToZero()
  {
    m_volumeMoments.setToZero();
    m_EBMoments.setToZero();
    for(int idir  = 0; idir < DIM; idir++)
    {
      m_normalPartialDeriv[idir].setToZero();
    }
    for(int iface  = 0; iface < 2*DIM; iface++)
    {
      m_faceMoments[iface].setToZero();
    }
  }
  ///
  void 
  IrregNode::
  setMomentsToRegular(const Real& a_dx)
  {
    //no eb for regular
    m_EBMoments.setToZero();
    for(int idir  = 0; idir < DIM; idir++)
    {
      m_normalPartialDeriv[idir].setToZero();
    }
    m_volumeMoments.setRegular(a_dx);
    for(int iface  = 0; iface < 2*DIM; iface++)
    {
      m_faceMoments[iface].setRegular(a_dx);
    }
  }
}
