
#ifndef _Proto_REALVECT_H_
#define _Proto_REALVECT_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>


#include "Proto_Point.H"
#include "Proto_IndexTM.H"
using std::vector;
namespace Proto
{
//template<typename T, int n>
//class IndexTM;

/// A double vector in DIM-dimensional space
/**
   The class RealVect is an implementation of a double vector in a
   DIM-dimensional space.
   RealVect values are accessed using the operator[] function, as for a normal
   C++ array.  In addition, the basic arithmetic operators have been overloaded
   to implement scaling and translation operations.
*/

  class RealVect
  {
  public:

    /**
       \name Constructors and Accessors
    */
    /*@{*/

    ///
    /**
       Construct a RealVect whose components are zero.
    */
    RealVect ()
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] = 0;
      }
    }

    ///
    /**
       The copy constructor.
    */
    RealVect (const RealVect& rhs)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] = rhs.vect[idir];
      }
    }

    ///
    /**
       Construct a RealVect from an IndexTM<double,DIM>.
    */
    RealVect(const IndexTM<double,DIM>& a_tm);

    ///
    /**
       Construct a RealVect from an Point by coercing each component
       from <tt>int</tt> to double.
    */
    RealVect (const Point & iv)
    {
      for (int d=0 ; d<DIM ; ++d)
      {
        vect[d] = (double)iv[d];
      }
    }

    ///
    /**
       The assignment operator.
    */
    RealVect& operator= (const RealVect& rhs)
    {
      for (int d=0 ; d<DIM ; ++d)
      {
        vect[d] = rhs.vect[d];
      }
    }

    ///
    /**
       Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
       RealVect.
    */
    inline
    double& operator[] (int i)
    {
      return vect[i];
    }

    ///
    /**
       Returns the <i>i</i>'th coordinate of the RealVect.
    */
    inline
    const double& operator[] (int i) const
    {
      return vect[i];
    }

    /**
       \name Arithmetic Operators
    */
    /*@{*/

    ///
    /**
       Modifies this RealVect by addition of a scalar to each component.
    */
    RealVect& operator+= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] += s;
      }
      return *this;
    }

    ///
    /**
       Returns a RealVect that is this RealVect with a scalar s added
       to each component.
    */
    RealVect operator+ (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] += s;
      }
      return retval;
    }

    ///
    /**
       Modifies this RealVect by component-wise addition by argument.
    */
    RealVect& operator+= (const RealVect& p)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] += p.vect[idir];
      }
      return *this;
    }

    ///
    /**
       Modifies this RealVect by subtraction of a scalar from each
       component.
    */
    RealVect& operator-= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] -= s;
      }
      return *this;
    }


    ///
    /**
       Modifies this RealVect by component-wise subtraction by argument.
    */
    RealVect& operator-= (const RealVect& p)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] -= p.vect[idir];
      }
      return *this;
    }

    ///
    /**
       Returns a RealVect that is this RealVect with a scalar s subtracted
       from each component.
    */
    RealVect operator- (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] -= s;
      }
      return retval;
    }

    ///
    /**
       Modifies this RealVect by multiplying each component by a scalar.
    */
    RealVect& operator*= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] *= s;
      }
      return *this;
    }


    ///
    /**

     */
    double dotProduct(const RealVect& a_rhs) const
    {
      double retval = 0;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval += vect[idir]*a_rhs.vect.[idir];
      }
      return retval;
    }


    ///
    /**
       Modifies this RealVect by component-wise multiplication by argument.
    */
    RealVect& operator*= (const RealVect& p)
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] *= p[idir];
      }
      return retval;
    }

//XXX  ///
//XXX  /**
//XXX     Returns component-wise product of this RealVect and argument.
//XXX  */
//XXX  RealVect operator* (const RealVect& p) const;

    ///
    /**
       Returns a RealVect that is this RealVect with each component
       multiplied by a scalar.
    */
    RealVect operator* (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] *= s;
      }
      return retval;
    }

    ///
    /**
       Modifies this RealVect by dividing each component by a scalar.
    */
    RealVect& operator/= (double s)
    {
      for(int idir = 0; idir < DIM; idir++)
      {
        vect[idir] /= s;
      }
      return *this;
    }


    ///
    /**
       Returns a RealVect that is this RealVect with each component
       divided by a scalar.
    */
    RealVect operator/ (double s) const
    {
      RealVect retval = *this;
      for(int idir = 0; idir < DIM; idir++)
      {
        retval.vect[idir] /= s;
      }
      return retval;
    }


    /*@}*/

    /**
       \name Other arithmetic operators
    */
    /*@{*/

    ///
    /**
       Modifies this RealVect by taking component-wise min with RealVect
       argument.
    */
    RealVect& min (const RealVect& p);

    ///
    /**
       Returns the RealVect that is the component-wise minimum of two
       argument RealVects.
    */
    friend inline RealVect min (const RealVect& p1,
                                const RealVect& p2);

    ///
    /**
       Modifies this RealVect by taking component-wise max with RealVect
       argument.
    */
    RealVect& max (const RealVect& p);

    ///
    /**
       Returns the RealVect that is the component-wise maximum of two
       argument RealVects.
    */
    friend inline RealVect max (const RealVect& p1,
                                const RealVect& p2);

    /*@}*/

    /**
       \name Unary operators
    */
    /*@{*/

    ///
    /**
       Unary plus -- for completeness.
    */
    RealVect operator+ () const;

    ///
    /**
       Unary minus -- negates all components of this RealVect.
    */
    RealVect operator- () const;

    ///
    /**
       Sum of all components of this RealVect.
    */
    double sum () const;

    ///
    /**
       sqrt(sum squares)
    */
    double vectorLength() const;

    ///
    /**
       sum squares--no square root
    */
    double radSquared() const;


    ///
    /**
       Component with the minimum value of this RealVect (returns 0 if they are all the same).
       a_doAbs : if true then take the absolute value before comparing
    */
    int minDir(const bool& a_doAbs) const
    {
      int retval = 0;
      if(a_doAbs)
      {
        double compval = std::abs(vect[0]);
        for(int idir = 1; idir < DIM; idir++)
        {
          if(std::abs(vect[idir]) < compval)
          {
            retval = idir;
            compval = std::abs(vect[idir]);
          }
        }
      }
      else
      {
        double compval = vect[0];
        for(int idir = 1; idir < DIM; idir++)
        {
          if(vect[idir] < compval)
          {
            retval = idir;
            compval = vect[idir];
          }
        }
      }
      return retval;
    }

    ///
    /**
       Component with the maximum value of this RealVect (returns 0 if they are all the same).
       a_doAbs : if true then take the absolute value before comparing
    */
    int maxDir(const bool& a_doAbs) const
    {
      int retval = 0;
      if(a_doAbs)
      {
        double compval = std::abs(vect[0]);
        for(int idir = 1; idir < DIM; idir++)
        {
          if(std::abs(vect[idir]) > compval)
          {
            retval = idir;
            compval = std::abs(vect[idir]);
          }
        }
      }
      else
      {
        double compval = vect[0];
        for(int idir = 1; idir < DIM; idir++)
        {
          if(vect[idir] > compval)
          {
            retval = idir;
            compval = vect[idir];
          }
        }
      }
      return retval;
    }


    /*@}*/


    /**
       \name Constants
    */
    /*@{*/

    ///
    /**
       Returns a basis vector in the given coordinate direction.<br>
       In 2-D:<br>
       BASISREALV(0) == (1.,0.);
       BASISREALV(1) == (0.,1.).<br>
       In 3-D:<br>
       BASISREALV(0) == (1.,0.,0.);
       BASISREALV(1) == (0.,1.,0.);
       BASISREALV(2) == (0.,0.,1.).<br>
       Note that the coordinate directions are based at zero.
    */
    friend RealVect BASISREALV(int dir);

    ///
    /**
       This is a RealVect all of whose components are equal to zero.
    */
    static const RealVect Zero = RealVect(0);

    ///
    /**
       This is a RealVect all of whose components are equal to one.
    */
    static const RealVect Unit= RealVect(1);

    /*@}*/

    /**
       \name Arithmetic friend functions
    */
    /*@{*/

    ///
    /**
       Returns a RealVect that is a RealVect <i>p</i> with
       a scalar <i>s</i> added to each component.
    */
    friend RealVect operator+ (double            s,
                               const RealVect& p);

    ///
    /**
       Returns <i>s - p</i>.
    */
    friend RealVect operator- (double            s,
                               const RealVect& p);

    ///
    /**
       Returns a RealVect that is a RealVect <i>p</i> with each component
       multiplied by a scalar <i>s</i>.
    */
    friend RealVect operator* (double            s,
                               const RealVect& p);
    ///
    /**
       Returns a RealVect that is a RealVect <i>p</i> with each component
       divided by a scalar <i>s</i>.
    */
    friend RealVect operator/ (double            s,
                               const RealVect& p);

    ///
    /**
       Returns component-wise sum of RealVects <i>s</i> and <i>p</i>.
    */
    friend RealVect operator+ (const RealVect& s,
                               const RealVect& p);

    ///
    /**
       Returns <i>s - p</i>.
    */
    friend RealVect operator- (const RealVect& s,
                               const RealVect& p);

    ///
    /**
       Returns component-wise product of <i>s</i> and <i>p</i>.
    */
    friend RealVect operator* (const RealVect& s,
                               const RealVect& p);
    ///
    /**
       Returns component-wise quotient <i>p / s</i>.
    */
    friend RealVect operator/ (const RealVect& s,
                               const RealVect& p);


    /*@}*/


  protected:

    /**
       The individual components of this RealVect.
    */
    double vect[DIM];

  };


  inline double& RealVect::operator[] (int i)
  {
    PROTO_ASSERT(i>=0 && i < DIM);
    return vect[i];
  }

  inline const double& RealVect::operator[] (int i) const
  {
    PROTO_ASSERT(i>=0 && i < DIM);
    return vect[i];
  }

  inline RealVect::RealVect (const RealVect &iv)
  {
    for(int idir = 0; idir < DIM; idir++)
    {
      vect[idir] = iv[idir];
    }
  }

  inline
  RealVect&
  RealVect::operator-= (double s)
  {
    for(int idir = 0; idir < DIM; idir++)
    {
      vect[idir] -= s;
    }
    return *this;
  }

  inline
  RealVect&
  RealVect::operator-= (const RealVect& p)
  {
    for(int idir = 0; idir < DIM; idir++)
    {
      vect[idir] -= p[idir];
    }
    return *this;
  }

  inline
  RealVect
  RealVect::operator+ () const
  {
    return RealVect(*this);
  }

  inline
  RealVect
  RealVect::operator- () const
  {
    RealVect retval = *this;
    retval *= -1.0;
    return retval;
  }

  inline
  double
  RealVect::sum () const
  {
    double retval = 0;
    for(int idir = 0; idir < DIM; idir++)
    {
      retval += vect[idir];
    }

    return retval;
  }

  inline
  double
  RealVect::vectorLength () const
  {
    double len = this->radSquared();
    len = sqrt(len);

    return len;
  }

  inline
  double
  RealVect::radSquared() const
  {
    double len = 0;
    for (int idir = 0; idir < DIM; idir++)
    {
      len = len + vect[idir]*vect[idir];
    }

    return len;
  }



  inline RealVect BASISREALV(int idir)
  {
    RealVect retval = RealVect::Zero;
    retval[idir] = 1;
    return rettval;
  }


}

#include "Proto_RealVectImplem.H"
#endif
