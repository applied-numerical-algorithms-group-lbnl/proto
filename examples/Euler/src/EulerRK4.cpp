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
    Bx dbx = dbx0.grow(NGHOST);
    BoxData<double,NUMCOMPS> U_ave(dbx);
    a_State.m_U.copyTo(U_ave);
    U_ave += a_DX.m_DU;
    for (auto iter = dbx.begin(); iter != dbx.end(); ++iter)
    {
        Point pt = *iter;
        if (!dbx0.contains(pt))
        {
            Point pt0 = pt;
            for (int dir = 0; dir < DIM; dir++)
            {
                pt0[dir] = dbx0.low()[dir] 
                    + (pt[dir] - dbx0.low()[dir])%dbx0.size(dir);
            }
            for (int k = 0; k < NUMCOMPS; k++)
            {
                U_ave(pt,k) = U_ave(pt0,k);
            }
        }
    }
    double velmax = m_eop(a_DX.m_DU,U_ave,a_State.m_dbx0);
    a_State.m_velSave = std::max(a_State.m_velSave,velmax);
    a_DX.m_DU *= a_dt;
  
    s_count += 1;
};
