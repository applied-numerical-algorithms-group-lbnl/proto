#include "EulerRK4.H"

int EulerRK4Op::s_count = 0;

void outToFile(const BoxData<double,NUMCOMPS>& a_U,const char* a_str);

EulerState::EulerState(const Bx& a_box)
{
    m_dbx0 = a_box;
    m_dbx0.print();
    m_U.define(m_dbx0.grow(NGHOST));
};

void 
EulerState::increment(const EulerDX& a_DX)
{m_U += a_DX.m_DU;};

void 
EulerDX::increment(double& a_weight,const EulerDX& a_incr)
{
    BoxData<double,NUMCOMPS> temp(a_incr.m_DU.box());
    (a_incr.m_DU).copyTo(temp);
    temp *= a_weight;
    m_DU += temp; 
};

void 
EulerDX::init(const EulerState& a_State)
{
    m_DU.define(a_State.m_dbx0);
    m_DU.setVal(0.);
};

void 
EulerDX::operator*=(const double& a_weight)
{m_DU *= a_weight;};

void 
EulerRK4Op::operator()(
                       EulerDX& a_DX,
                       double a_time,
                       double a_dt,
                       EulerState& a_State)
{
    Bx dbx0 = a_State.m_dbx0;
    Bx domain = dbx0;
    Bx dbx = dbx0.grow(NGHOST);
    BoxData<double,NUMCOMPS> U_ave(dbx);
    a_State.m_U.copyTo(U_ave);
    U_ave += a_DX.m_DU;

    int numghost = NGHOST;
    int numcomps = NUMCOMPS;
    for(int idir = 0; idir < DIM; idir++)
    {
      Point dirlo = -Point::Basis(idir);
      Point dirhi =  Point::Basis(idir);
      Bx edgebxlo = domain.edge(dirlo);
      Bx edgebxhi = domain.edge(dirhi);

      Bx dstbxlo = edgebxlo.extrude(idir,-(numghost-1));
      Bx dstbxhi = edgebxhi.extrude(idir, (numghost-1));

      //these are swapped because you have to move in different 
      //directions to get the periodic image
      Point shifthi = dirlo*domain.size(idir);
      Point shiftlo = dirhi*domain.size(idir);

      Bx srcbxlo = dstbxlo.shift(shiftlo);
      Bx srcbxhi = dstbxhi.shift(shifthi);

      U_ave.copy(U_ave, srcbxlo, 0, dstbxlo, 0, numcomps);
      U_ave.copy(U_ave, srcbxhi, 0, dstbxhi, 0, numcomps);

    }


    double velmax = EulerOp::step(a_DX.m_DU,U_ave,a_State.m_dbx0);
    a_State.m_velSave = std::max(a_State.m_velSave,velmax);
    a_DX.m_DU *= a_dt;
  
    s_count += 1;
};
