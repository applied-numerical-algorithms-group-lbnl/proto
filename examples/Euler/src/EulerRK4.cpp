#include "EulerRK4.H"

int EulerRK4Op::s_count = 0;

void outToFile(const BoxData<double,NUMCOMPS>& a_U,const char* a_str);

EulerState::EulerState(const Box& a_box)
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
    Box dbx0 = a_State.m_dbx0;
    Box domain = dbx0;
    Box dbx = dbx0.grow(NGHOST);
    BoxData<double,NUMCOMPS> U_ave(dbx);
    a_State.m_U.copyTo(U_ave);
    U_ave += a_DX.m_DU;

//    int numghost = NGHOST;
//    int numcomps = NUMCOMPS;
//    for(int idir = 0; idir < DIM; idir++)
//    {
//      Point dirlo = -Point::Basis(idir);
//      Point dirhi =  Point::Basis(idir);
//      Box edgebxlo = domain.edge(dirlo);
//      Box edgebxhi = domain.edge(dirhi);
//
//      Box dstbxlo = edgebxlo.extrude(idir,-(numghost-1));
//      Box dstbxhi = edgebxhi.extrude(idir, (numghost-1));
//
//      //these are swapped because you have to move in different 
//      //directions to get the periodic image
//      Point shifthi = dirlo*domain.size(idir);
//      Point shiftlo = dirhi*domain.size(idir);
//
//      Box srcbxlo = dstbxlo.shift(shiftlo);
//      Box srcbxhi = dstbxhi.shift(shifthi);
//
//      U_ave.copy(U_ave, srcbxlo, 0, dstbxlo, 0, numcomps);
//      U_ave.copy(U_ave, srcbxhi, 0, dstbxhi, 0, numcomps);
//
//    }


    Reduction<double>& rxn = a_State.m_Rxn;
    rxn.reset();
    EulerOp::step(a_DX.m_DU,U_ave,a_State.m_dbx0, rxn, false, true);
    a_State.m_velSave = std::max(a_State.m_velSave, rxn.fetch());
    a_DX.m_DU *= a_dt;
  
    s_count += 1;
};
