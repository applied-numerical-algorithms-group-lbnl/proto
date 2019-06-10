
#include "Proto.H"
#include "Proto_LevelData.H"
#include "EulerOp.H"

  
int main(int argc, char* argv[])
{
  Box domain(-16*Point::Unit(), 47*Point::Unit());

  std::array<bool, DIM> periodic;
  for(int idir = 0; idir < DIM; idir++) periodic[idir]=true;

  DisjointBoxLayout   dbl(domain, 16, periodic);

  LevelData<BoxData<double, NUMCOMPS>> U(dbl, NGHOST*Point::Unit());
  LevelData<BoxData<double, NUMCOMPS>> RHS(dbl, Point::Zero());
  
  for(unsigned int i=0; i<dbl.size(); i++)
    {
      auto u = U[i];
      auto rhs = RHS[i];
      Box rbox = dbl[i];
      double wave = EulerOp::step(rhs, u, rbox);
    }

    for(unsigned int i=0; i<dbl.size(); i++)
    {
      auto u = U[i];
      auto rhs = RHS[i];
      Box rbox = dbl[i];
      double wave = EulerOp::step(rhs, u, rbox);
    }

    return 0;
}
  
      
