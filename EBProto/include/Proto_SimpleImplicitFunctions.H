
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
    virtual bool entireBoxRegular(const Box& a_box,
                                  const double& a_dx)
    {
      vector<RealVect> cornerlocs = getCornerLocations(a_box, a_dx);
      bool allNeg = true;

      IndexTM<int, DIM> deriv;
      deriv.setAll(0);
      for(int icorn = 0; icorn < cornerlocs.size(); icorn++)
      {
        IndexTM<double, DIM> itmloc = getIndexTM(cornerlocs[icorn]);
        double value = value(deriv, itmloc);
        if(value > 0) 
        {
          allNeg = false;
        }
      }
      return allNeg;
    }

    //overwrite if you have a fast way of knowing if a whole region is regular
    virtual bool entireBoxCovered(const Box& a_box,
                                  const double& a_dx)
    {
      vector<RealVect> cornerlocs = getCornerLocations(a_box, a_dx);
      bool allPos = true;

      IndexTM<int, DIM> deriv;
      deriv.setAll(0);
      for(int icorn = 0; icorn < cornerlocs.size(); icorn++)
      {
        IndexTM<double, DIM> itmloc = getIndexTM(cornerlocs[icorn]);
        double value = value(deriv, itmloc);
        if(value < 0) 
        {
          allPos = false;
        }
      }
      return allPos;
    }

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
      Point hipt = a_box.high + Point::Ones(1);
#if DIM==2
      int numcorners = 4;
      cornerlocs.resize(numcorners);

      Point lohi = lopt;
      lohi[1] = hipt[1];   
      Point hilo = hipt;
      hilo[1] = lopt[1];   

      cornerlocs[0] = a_dx*lopt;
      cornerlocs[1] = a_dx*hipt;
      cornerlocs[2] = a_dx*lohi;
      cornerlocs[3] = a_dx*hilo;
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
      hilolo[0] = hipt[0];   


      Point lohilo = lopt;
      lolohi[1] = hipt[1];   
      Point hilohi = hipt;
      hihilo[1] = lopt[1];   

      cornerlocs[0] = a_dx*lopt;
      cornerlocs[1] = a_dx*hipt;

      cornerlocs[2] = a_dx*lolohi;
      cornerlocs[3] = a_dx*hihilo;

      cornerlocs[4] = a_dx*hilolo;
      cornerlocs[5] = a_dx*lohihi;

      cornerlocs[6] = a_dx*hilohi;
      cornerlocs[7] = a_dx*lohilo;
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

      for(int idir = 0; idir < SpaceDim; idir++)
      {
        CH_assert(a_ABE[idir] > 1.e-10);
        m_A2[idir] = a_ABE[idir]*a_ABE[idir];
      }
    }

    /// Destructor
    virtual ~SimpleEllipsoidIF()
    {
    }

    double
    value(const IndexTM<Real,SpaceDim>& a_point) const
    {
  
      RealVect X;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      Real retval = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        retval += (X[idir]*X[idir])/m_A2[idir];
      }

      retval -= m_R*m_R;
      if(!m_inside)
        retval = -retval;
      return retval;
    }

    virtual double value(const IndexTM<int,SpaceDim> & a_partialDerivative,
                         const IndexTM<double,SpaceDim>& a_point) const
    {
      int order = a_partialDerivative.sum();
      double retval = LARGEREALVAL;

      RealVect X;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      if (order == 0)
      {
        retval = value(a_point);
      }
      else if (order == 1)
      {
        bool found = false;
        for(int idir = 0; idir < SpaceDim; idir++)
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
        for(int idir = 0; idir < SpaceDim; idir++)
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
  class SimpleSphere: public SimpleGeometry
  {
  public:
    ///
    SimpleSphereIF(const RealVect& a_X0,
                   const double  & a_R,
                   const bool    & a_inside)
    {
      m_X0=       a_X0;    
      m_R=        a_R;     
      m_inside=   a_inside;

        m_A2[idir] = a_ABE[idir]*a_ABE[idir];
      }
    }

    /// Destructor
    virtual ~SimpleSphereIF()
    {
    }

    double
    value(const IndexTM<Real,SpaceDim>& a_point) const
    {
  
      RealVect X;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      Real retval = 0;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        retval += (X[idir]*X[idir]);
      }

      retval -= m_R*m_R;
      if(!m_inside)
        retval = -retval;
      return retval;
    }

    virtual double value(const IndexTM<int,SpaceDim> & a_partialDerivative,
                         const IndexTM<double,SpaceDim>& a_point) const
    {
      int order = a_partialDerivative.sum();
      double retval = LARGEREALVAL;

      RealVect X;
      for(int idir = 0; idir < SpaceDim; idir++)
      {
        X[idir] = a_point[idir] - m_X0[idir];
      }

      if (order == 0)
      {
        retval = value(a_point);
      }
      else if (order == 1)
      {
        bool found = false;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_partialDerivative[idir] == 1)
          {
            found = true;
            retval = 2*X[idir];
          }
        }
       }
      else if (order == 2)
      {
      
        bool found = false;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
          if(a_partialDerivative[idir] == 2)
          {
            found = true;
            retval = 2;
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
      SimpleSphereIF* SimpleSpherePtr = new SimpleSphereIF(m_ABE,m_X0, m_R, m_inside);
      return static_cast<BaseIF*>(SimpleSpherePtr);
    }

  protected:
    RealVect m_ABE;
    RealVect m_X0;
    double   m_R;
    bool     m_inside;
    RealVect m_A2;

  private:
    SimpleSphereIF();

  };

}
#endif
