#include "Proto.H"
#include <map>

#define NUMCOMPS 1
#define NGHOST 2



namespace Proto {
class BoxDataDoubleDecay
{
public:
  const Box m_box;
  const double* m_data;
  unsigned int m_C;
  template <unsigned int C>
  BoxDataDoubleDecay(BoxData<double, C>& m_arg):
    m_box(m_arg.box()),
    m_data(m_arg.data()),
    m_C(C) {}
};

template <typename ...Args>
int PROTOX_START()
{
  //last case, wrap up
  return 0;
}

template <typename ...Args>  
int PROTOX_START(BoxDataDoubleDecay First, Args&&...  m_rest)
{

  // handle BoxData<double, C> case

  return 1+PROTOX_START(std::forward<Args>(m_rest)...);
};

template <typename ...Args>
int PROTOX_START(double First, Args&&...  m_rest)
{

  // handle floating point case

  return 1+PROTOX_START(std::forward<Args>(m_rest)...);
};

template <typename ...Args>  
int PROTOX_START(int First, Args&&...  m_rest)
{

  // handle integer case

  return 1+PROTOX_START(std::forward<Args>(m_rest)...);
};

template <typename ...Args> 
int PROTOX_START(Box First, Args&&...  m_rest)
{

  // handle Box case

  return 1+PROTOX_START(std::forward<Args>(m_rest)...);
};
};


using namespace Proto;

int main(int argc, char *argv[])
{
  
  int nx = 32;
  Box bx(Point::Ones(-NGHOST),Point::Ones(nx+NGHOST));
  Box bx0(Point::Ones(1),Point::Ones(nx-1));
  BoxData<double,NUMCOMPS> u(bx);
  BoxData<double,NUMCOMPS> Lofu(bx0);

  Box applyBox(Point::Ones(3), Point::Ones(nx-3));
  Array<double,3> error;
  double gamma=1.4;

  int argCount = PROTOX_START(u, Lofu, 17, gamma, bx, bx0);
 
  std::cout<<"arcgCount="<<argCount<<std::endl;

  return argCount;
}
         
                        
