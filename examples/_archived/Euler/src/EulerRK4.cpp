#include "EulerRK4.H"

int EulerRK4Op::s_count = 0;

EulerState::EulerState(const Box& a_box,
                       const double a_dx,
                       const double a_gamma):
    m_dx(a_dx),
    m_gamma(a_gamma)
{
    m_dbx0 = a_box;
    Box dbx = m_dbx0.grow(NGHOST);
    m_U.define(m_dbx0);
    Point ptValid = Point::Zeros();
    // Setup for periodic BCs. Good for single Box data only.
    for (int dir = 0; dir < DIM; dir++)
      {
       int size = m_dbx0.size(dir);
       Box bValid = m_dbx0.grow(ptValid*NGHOST);
       m_shift[2*dir]=Point::Basis(dir)*size;
       m_bdry[2*dir] = bValid&dbx.shift(-m_shift[2*dir]);
       m_shift[2*dir+1] = Point::Basis(dir)*(-size);
       m_bdry[2*dir+1] = bValid&dbx.shift(-m_shift[2*dir+1]);
       ptValid = ptValid + Point::Basis(dir);
      }
};
void 
EulerState::increment(const EulerDX& a_DX)
{
  m_U += a_DX.m_DU;
};
void 
EulerDX::increment(double& a_weight,const EulerDX& a_incr)
{
  BoxData<double,NUMCOMPS> temp(a_incr.m_DU.box());
  (a_incr.m_DU).copyTo(temp);
  temp *= a_weight;
  m_DU += temp;
};

void 
EulerDX::init(EulerState& a_State)
{
    m_DU.define(a_State.m_dbx0);
    m_DU.setVal(0.);
    m_box = a_State.m_dbx0;
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
    (a_State.m_U).copyTo(U_ave);
    
     
    U_ave += a_DX.m_DU;
    for(int dir = 0; dir < DIM; dir++)
      {
        //std::cout << 2*dir << " , " << a_State.m_bdry[2*dir] << " , "<< a_State.m_shift[2*dir] << std::endl;
        //std::cout << U_ave.box() << " , " << U_ave.box().shift(a_State.m_shift[0]) << std::endl;
        U_ave.copyTo(U_ave,a_State.m_bdry[2*dir], a_State.m_shift[2*dir]);
        //std::cout << 2*dir+1 << std::endl;
        U_ave.copyTo(U_ave,a_State.m_bdry[2*dir+1], a_State.m_shift[2*dir+1]);        
      }
    Reduction<double>& rxn = a_State.m_Rxn;
    EulerOp::step(a_DX.m_DU,U_ave,a_State.m_dbx0, a_State.m_dx, a_State.m_gamma, rxn, true, true);
    a_DX*=a_dt;
    s_count += 1;
};
