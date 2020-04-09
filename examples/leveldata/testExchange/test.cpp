#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <memory>
#include <stdio.h>
#include <fstream>
#include "Proto.H"
#include <iostream>
#include <cstring>
#include <memory>
#define MAXBOXSIZE 64
using namespace std;
using namespace Proto; 
int main(int argc, char* argv[])
{
  int domainSize = 64;
  Box domain(Point::Zeros(),Point::Ones()*(domainSize -1));

  array<bool,DIM> per;
  for(int idir = 0; idir < DIM; idir++) per[idir]=true;

  DisjointBoxLayout bl(domain,min(MAXBOXSIZE,domainSize),per);

  LevelData<BoxData<double> > phi(bl,Point::Ones());
  for (int i = 0; i < bl.size();i++)
    {
      BoxData<double>& phiPatch = phi[i];
      phiPatch.setVal(1.);
    }
  phi.exchange();
}
