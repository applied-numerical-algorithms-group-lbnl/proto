
#ifndef ___Proto_SimpleImplicitFunctions__H__
#define ___Proto_SimpleImplicitFunctions__H__

#include "Proto_RealVect.H"
#include "Proto_IndexTM.H"
#include "Proto_BaseIF.H"
namespace Proto
{
  ///truly simple geometries can have some simplified bits about what is all covered or regular
  class SimpleGeometry: public BaseIF
  {
  public:


    IndexTM<double,DIM> getIndexTM(const RealVect& a_rv)
    {
      IndexTM<double, DIM> retval;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval[idir] = a_rv[idir];
      }
      return retval;
    }

    vector<RealVect> getCornerLocations(const Box    & a_box,
                                        const double & a_dx)
    {
      vector<RealVect> cornerlocs;
      Point lopt = a_box.low();
      Point hipt = a_box.high() + Point::Ones(1);
#if DIM==2
      int numcorners = 4;
      cornerlocs.resize(numcorners);

      Point lohi = lopt;
      lohi[1] = hipt[1];   
      Point hilo = hipt;
      hilo[1] = lopt[1];   

      for(int idir = 0; idir < DIM; idir++)
      {
        cornerlocs[0][idir] = a_dx*lopt[idir];
        cornerlocs[1][idir] = a_dx*hipt[idir];
        cornerlocs[2][idir] = a_dx*lohi[idir];
        cornerlocs[3][idir] = a_dx*hilo[idir];
      }
#else
      int numcorners = 8;

      cornerlocs.resize(numcorners);
      Point lolohi = lopt;
      lolohi[2] = hipt[2];   
      Point hihilo = hipt;
      hihilo[2] = lopt[2];   

      Point lohihi = hipt;
      lohihi[0] = lopt[0];   
      Point hilolo = lopt;



      Point lohilo = lopt;
      lolohi[1] = hipt[1];   
      Point hilohi = hipt;
      hihilo[1] = lopt[1];   

      for(int idir = 0; idir < DIM; idir++)
      {
        cornerlocs[0][idir] =   a_dx*lopt[idir];
        cornerlocs[1][idir] =   a_dx*hipt[idir];
        cornerlocs[2][idir] = a_dx*lolohi[idir];
        cornerlocs[3][idir] = a_dx*hihilo[idir];
        cornerlocs[4][idir] = a_dx*hilolo[idir];
        cornerlocs[5][idir] = a_dx*lohihi[idir];
        cornerlocs[6][idir] = a_dx*hilohi[idir];
        cornerlocs[7][idir] = a_dx*lohilo[idir];
      }
#endif


      return cornerlocs;
    }
  };

///
/**
   ((x-xo)/A)^2 + ((y-y0)/B)^2 + ((z-z0)/C)^2 - r^2 = f(x,y,z)
*/
  class SimpleEllipsoidIF: public SimpleGeometry
  {
  public:
    ///
    SimpleEllipsoidIF(const RealVect& a_ABE,
                      const RealVect& a_X0,
                      const double  & a_R,
                      const bool    & a_inside)
    {
      m_ABE=      a_ABE;   
      m_X0=       a_X0;    
      m_R=        a_R;     
      m_inside=   a_inside;

      for(int idir = 0; idir < DIM; idir++)
      {
        PR_assert(a_ABE[idir] > 1.e-10);
        m_A2[idir] = a_ABE[idir]*a_ABE[idir];
      }
    }

    bool sameQuadrant(vector<RealVect> a_cornerlocs)
    {
      vector<int> quadrant(a_cornerlocs.size());
      for(int iloc = 0; iloc < a_cornerlocs.size(); iloc++)
      {
        const RealVect& loc = a_cornerlocs[iloc];
        quadrant[iloc] = -1;
#if DIM==2
        if((loc[0] > m_X0[0]) && (loc[1] > m_X0[1]))
        {
          quadrant[iloc] = 0;
        }
        if((loc[0] > m_X0[0]) && (loc[1] < m_X0[1]))
        {
          quadrant[iloc] = 1;
        }
        if((loc[0] < m_X0[0]) && (loc[1] < m_X0[1]))
        {
          quadrant[iloc] = 2;
        }
        if((loc[0] < m_X0[0]) && (loc[1] > m_X0[1]))
        {
          quadrant[iloc] = 3;
        }
#else
        if((loc[0] > m_X0[0]) && (loc[1] > m_X0[1]) && (loc[2] > m_X0[2]))
        {
          quadrant[iloc] = 0;
        }
        if((loc[0] > m_X0[0]) && (loc[1] < m_X0[1]) && (loc[2] > m_X0[2]))
        {
          quadrant[iloc] = 1;
        }
        if((loc[0] < m_X0[0]) && (loc[1] < m_X0[1]) && (loc[2] > m_X0[2]))
        {
          quadrant[iloc] = 2;
        }
        if((loc[0] < m_X0[0]) && (loc[1] > m_X0[1]) && (loc[2] > m_X0[2]))
        {
          quadrant[iloc] = 3;
        }
        if((loc[0] > m_X0[0]) && (loc[1] > m_X0[1]) && (loc[2] < m_X0[2]))
        {
          quadrant[iloc] = 4;
        }
        if((loc[0] > m_X0[0]) && (loc[1] < m_X0[1]) && (loc[2] < m_X0[2]))
        {
          quadrant[iloc] = 5;
        }
        if((loc[0] < m_X0[0]) && (loc[1] < m_X0[1]) && (loc[2] < m_X0[2]))
        {
          quadrant[iloc] = 6;
        }
        if((loc[0] < m_X0[0]) && (loc[1] > m_X0[1]) && (loc[2] < m_X0[2]))
        {
          quadrant[iloc] = 7;
        }
#endif
      }
      bool allSame = true;
      for(int iloc = 1; iloc < quadrant.size(); iloc++)
      {
        if(quadrant[iloc] != quadrant[0])
        {
          allSame = false;
        }
      }
      return allSame;
    }


    virtual bool entireBoxRegular(const Box& a_box,
                                  const double& a_dx)
    {
      vector<RealVect> cornerlocs = getCornerLocations(a_box, a_dx);
      bool allSameQuadrant = sameQuadrant(cornerlocs);
      bool retval = false;
      if(allSameQuadrant)
      {
        bool allNeg = true;

        IndexTM<int, DIM> deriv;
        deriv.setAll(0);
        for(int icorn = 0; icorn < cornerlocs.size(); icorn++)
        {
          IndexTM<double, DIM> itmloc = getIndexTM(cornerlocs[icorn]);
          double funcvalue = value(deriv, itmloc);
          if(funcvalue > 0) 
          {
            allNeg = false;
          }
          retval = allNeg;
        }
      }
      return retval;
    }
    //overwrite if you have a fast way of knowing if a whole region is regular
    virtual bool entireBoxCovered(const Box& a_box,
                                  const double& a_dx)
    {
      vector<RealVect> cornerlocs = getCornerLocations(a_box, a_dx);
      bool allSameQuadrant = sameQuadrant(cornerlocs);
      bool retval = false;
      if(allSameQuadrant)
      {

        bool allPos = true;

        IndexTM<int, DIM> deriv;
        deriv.setAll(0);
        for(int icorn = 0; icorn < cornerlocs.size(); icorn++)
        {
          IndexTM<double, DIM> itmloc = getIndexTM(cornerlocs[icorn]);
          double funcvalue = value(deriv, itmloc);
          if(funcvalue < 0) 
          {
            allPos = false;
          }
        }
        retval = allPos;
      }
      return retval;
    }
    /// Destructor
    virtual ~SimpleEllipsoidIF()
    {
    }

    double
    value(const IndexTM<double,DIM>& a_point) const
    {
  
      RealVect X;
      for(int idir = 0; idir < DIM; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      double retval = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval += (X[idir]*X[idir])/m_A2[idir];
      }

      retval -= m_R*m_R;
      if(!m_inside)
        retval = -retval;
      return retval;
    }

    virtual double value(const IndexTM<int,DIM> & a_partialDerivative,
                         const IndexTM<double,DIM>& a_point) const
    {
      int order = a_partialDerivative.sum();
      double retval = LARGEREALVAL;

      RealVect X;
      for(int idir = 0; idir < DIM; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      if (order == 0)
      {
        retval = value(a_point);
        return retval; //have to return here to avoid double neg
      }
      else if (order == 1)
      {
        bool found = false;
        for(int idir = 0; idir < DIM; idir++)
        {
          if(a_partialDerivative[idir] == 1)
          {
            found = true;
            retval = 2*X[idir]/m_A2[idir];
          }
        }
      }
      else if (order == 2)
      {
      
        bool found = false;
        for(int idir = 0; idir < DIM; idir++)
        {
          if(a_partialDerivative[idir] == 2)
          {
            found = true;
            retval = 2/m_A2[idir];
          }
        }
        if(!found) //mixed deriv
        {
          retval = 0;
        }
      }
      else
      {
        retval = 0.0;
      }

      if(!m_inside)
        retval = -retval;
      return retval;
    }


    virtual BaseIF* newImplicitFunction() const
    {
      SimpleEllipsoidIF* SimpleEllipsoidPtr = new SimpleEllipsoidIF(m_ABE,m_X0, m_R, m_inside);
      return static_cast<BaseIF*>(SimpleEllipsoidPtr);
    }

  protected:
    RealVect m_ABE;
    RealVect m_X0;
    double   m_R;
    bool     m_inside;
    RealVect m_A2;

  private:
    SimpleEllipsoidIF();

  };



///
/**
   (x-xo)^2 + (y-y0)^2 + (z-z0)^2 - r^2 = f(x,y,z)

*/
  class SimpleSphereIF: public SimpleEllipsoidIF
  {
  public:
    ///
    SimpleSphereIF(const RealVect& a_X0,
                   const double  & a_R,
                   const bool    & a_inside):
      SimpleEllipsoidIF(RealVect::Unit(), a_X0, a_R, a_inside)
    {
    }

    /// Destructor
    virtual ~SimpleSphereIF()
    {
    }


  private:
    SimpleSphereIF();

  };

}
#endif
