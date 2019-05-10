#ifndef _AGGSTENCIL_H_
#define _AGGSTENCIL_H_

#include "EBProto.H"

namespace Proto
{
  ///replacement for Chombo's vofstencil and facestencil
  template<CENTERING cent, typename data_t>
  class LocalStencil
  {
  public:
    struct entry_t
    {
      EBIndex<cent> m_vof;
      data_t        m_weight;
      entry_t(const EBIndex<cent>& a_vof,
              const data_t       & a_weight)
      {
        m_vof = a_vof;
        m_weight = a_weight;
      }
    };
    vector<element_t>  m_entries;


    ///
    LocalStencil()
    {}
    
    ///
    ~LocalStencil()
    {}
    
    ///
    inline void clear()
    {
      m_entries.resize(0);
    }
    
    ///
    int size() const
    {
      return m_entries.size();
    }
    
    ///this api should be stable
    inline void add(const entry_t& a_entry)
    {
      add(a_entry.m_vof, a_entry.m_weight);
    }

    ///this might change (if a variable number argument is added, for example)
    inline void add(const EBIndex<cent> a_vof, const data_t& a_weight)
    {
      bool alreadyhere = false;
      for (int ivof = 0; ivof < m_vofs.size(); ivof++)
      {
        if(m_entries[ivof].m_vof == a_vof)
        {
          alreadyhere = true;
          m_entries[ivof].m_weight += weight;
        }
      }
      if (!alreadyhere)
      {
        m_entries.push_back(entry_t(a_vof, a_weight));
      }
    }
   
    ///add full stencil
    LocalStencil& operator+=(const LocalStencil<cent>& a_input)
    {
      if(&a_input != this)
      {
        for(int ient = 0; ient < a_input.size(); ient++)
        {
          add(a_input.m_entries[ient]);
        }
      }
      return *this;
    }

    ///scaling operator
    void operator*=  (const data_t & a_scaling)
    {
      for(int ient = 0; ient < m_entries.size(); ient++)
      {
        m_entries[ient].m_weight *= a_scaling;
      }
    }
  };


  /// 
  /**
  */
  template <CENTERING srcCenter, CENTERING dstCenter, typename data_t>
  class AggStencil
  {
  public:
    ///
    /**
    */
    AggStencil(const vector< EBIndex<     dstCenter        > > & a_dstVoFs,
               const vector< LocalStencil<srcCenter, data_t> > & a_stencil,
               const EBGraph                                   & a_graph,
               const Box                                       & a_validBox,
               const Point                                     & a_srcGhost,
               const Point                                     & a_dstGhost);


    ///
    ~AggStencil()
    {
    }

     ///
    /**
       lphi = L(phi)
     */
    void apply(EBBoxData<      dstCenter, data_t, 1>     & a_lph,
               const EBBoxData<srcCenter, data_t, 1>     & a_phi,
               const data_t                              & a_scale);

    ///EBDataLoc defined in Proto_EBBoxData.H
    typedef vector<pair<EBDataLoc, data_t> > stencil_t;

  protected:

    vector<stencil_t>    m_ebstencil;
    vector<EBDataLoc>    m_dstAccess;

  private:
    /// disallowed operators.  
    AggStencil();
    AggStencil& operator=(const AggStencil& stenin);
    AggStencil(const AggStencil& stenin);

  };
}
#include "Proto_AggStencilImplem.H"

#endif
