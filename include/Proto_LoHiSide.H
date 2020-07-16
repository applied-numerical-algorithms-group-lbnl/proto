#ifndef _Proto_LOHISIDE_H_
#define _Proto_LOHISIDE_H_

namespace Proto
{
/// Encapsulation of high or low side.
/**
   Encapsulation of high or low side.
*/
  class Side
  {
  public:
    ///
    /**
       Lo is the low side.<br>
       Hi is the high side.
    */
    enum LoHiSide
    {
      Invalid=-1,
      Lo=0,
      Hi,
      NUMSIDES
    };

    /**
       Constructor.
    */
    Side()
    {
    }

    /**
       Destructor.
    */
    ~Side()
    {
    }

    ///
    inline static Side::LoHiSide flip(const Side::LoHiSide& a_side)
    {
      return (a_side == Side::Lo) ? Side::Hi : Side::Lo;
    }

  };

///
  inline int sign(const Side::LoHiSide& a_side)
  {
    return (a_side == Side::Lo) ? -1 : 1;
  }

///
  inline Side::LoHiSide flip(const Side::LoHiSide& a_side)
  {
    //PR_assert((a_side == Side::Lo) || (a_side == Side::Hi));

    return (a_side == Side::Lo) ? Side::Hi : Side::Lo;
  }
  
/// Iterator for low and high side.
/**
   SideIterator iterates over the values Side::Lo and Side::Hi.
   There are no adjustable parameters whatsoever.
*/
  class SideIterator
  {
  public:
    ///
    /**
       Default, and only, constructor.
    */
    SideIterator()
    {
      m_current = -1;
      reset();
    }

    ///
    /**
       Default destructor.
    */
    ~SideIterator()
    {}

    /**
       \name Modification functions
    */
    /*@{*/

    ///
    /**
       Set the iterator to the first side.
    */
    inline void begin()
    {
      m_current = 0;
    }

    ///
    /**
       Set the iterator to the first side.
    */
    inline void reset()
    {
      begin();
    }

    ///
    /**
       Advance the iterator to the next face.
    */
    inline void next()
    {
      ++m_current;
    }

    ///
    /**
       Advance the iterator to the next face.
    */
    inline void operator ++ ()
    {
      ++m_current;
    }

    /*@}*/

    /**
       \name Access functions
    */
    /*@{*/

    ///
    /**
       Returns the value of the side for the current iteration.
    */
    inline Side::LoHiSide operator () () const
    {
      switch (m_current)
      {
      case 0:
        return Side::Lo;
        //break;

      case 1:
        return Side::Hi;
        //break;

      default:
        return Side::Invalid;
        //break;
      }
    }


    ///
    /**
       Returns <tt>true</tt> if there is a next side, <tt>false</tt> if not.
    */
    inline bool ok() const
    {
      return ((m_current > -1) && (m_current < Side::NUMSIDES));
    }


  private:
    int m_current;
  };

}
#endif
