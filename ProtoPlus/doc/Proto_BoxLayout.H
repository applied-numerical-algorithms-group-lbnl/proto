#pragma once

namespace Proto {

// forward declarations
class DataIterator;
class LayoutIterator;

/// Disjoint Box Layout
/**
  Implicitly assumes disjointness and fixed Box size.
 */
class BoxLayout
{
public:
    inline BoxLayout(Box a_domain, //physical domain
            unsigned int a_patchSize = PROTO_PATCHSIZE, //fixed size of each patch (in all directions; can make this a point for anisotropy)
            Point a_periodic = Point::Ones())
    {
        Box bitbox = a_domain.coarsen(a_patchSize);
        BoxData<bool> bitmap(bitbox, true);
        // do load balancing here -> build patchRanks map
        define(patchRanks, bitmap, a_patchSize, a_periodic);
    }

    inline void define(std::map<Point, unsigned int> a_patchRanks,  
            BoxData<bool> a_bitmap,                      
            unsigned int a_patchSize = PROTO_PATCHSIZE, 
            Point a_periodic = Point::Ones());               
    {
        //m_rank obtained from MPI
        m_patchSize = a_patchSize;
        m_periodic = a_periodic;
        m_patchRanks = a_patchRanks;
        m_bitmap = a_bitmap;
    }

    inline DataIterator dataIterator() const;
    inline LayoutIterator layoutIterator() const;
    inline Box operator[](const DataIterator& a_iter) const;
    inline Box operator[](const LayoutIterator& a_iter) const;

    // coarsen and refine are really easy; the only data member that changes is m_patchSize
    inline BoxLayout coarsen(unsigned int a_ratio) const;
    inline BoxLayout refine(unsigned int a_ratio) const;

private:
    unsigned int m_rank;             ///< This parallel process rank
    unsigned int m_patchSize;        ///< Size of all Boxes in the array. 
    // This information might need to live in LevelData; I don't think the separate ProblemDomain object used in Chombo is necessary.
    Point m_periodic;                ///< Defines periodicity in each direction. m_periodic[dir] == implies non-periodic, else the domain is periodic.
    std::map<Point, unsigned int> m_patchRanks; ///< Maps points in m_bitmap to a rank. if m_patchRanks[point] == m_rank -> Box(point,point).refine(m_patchSize) lives here.
    BoxData<bool> m_bitmap;         ///< Each point represents a possible patch in the domain. Domain box is m_domain.coarsen(m_patchSize).
}

/// Local Data Iterator
/**
  Used to iterate over all the data on this proc; ignores data on other processors. Can be used by any LevelData defined with the same BoxLayout (though possibly different patchSizes)
 */
class DataIterator
{
public:
    inline LayoutIterator(BoxData<bool>& a_bitmap,
            std::map<Point, int>& a_patchRanks,
            unsigned int a_rank)
    {
        m_bitmap = a_bitmap;
        m_patchRanks = a_patchRanks;
        m_rank = a_rank;
        this->begin();
    }
    inline Point operator*()
    {
        return *m_iter;
    }
    inline void begin()
    {
        m_iter = a_bitmap.box().begin();
        while ((!m_bitmap(*m_iter)) || (m_patchRanks[*m_iter] != m_rank))
        {
            m_iter++;
            if (m_iter.done()){break;}
        }
    }
    inline void operator++()
    {
        m_iter++;
        while ((!m_bitmap(*m_iter)) || (m_patchRanks[*m_iter] != m_rank))
        {
            m_iter++;
            if (m_iter.done()){break;}
        }
    }
    inline bool done(){m_iter.done();}
private:
    const BoxData<bool>& m_bitmap;
    const std::map<Point, unsigned int>& m_patchRanks;
    BoxIterator m_iter;
    unsigned int m_rank;
}
/// Layout Iterator
/**
  Iterates over the entire BoxLayout, making no distinction vis a vis local / remove processes
 */
class LayoutIterator
{
    public:
        inline LayoutIterator(BoxData<bool>& a_bitmap)
        {
            m_bitmap = a_bitmap;
            begin();
        }
        inline Point operator*()
        {
            return *m_iter;
        }
        inline void begin()
        {
            m_iter = a_bitmap.box().begin();
            while (!m_bitmap(*m_iter)){m_iter++;}
        }
        inline void operator++()
        {
            m_iter++;
            while (!m_bitmap(*m_iter)){m_iter++;}
        }
        inline bool done(){m_iter.done();}
    private:
        const BoxData<bool>& m_bitmap;
        BoxIterator m_iter;
}

} //end Proto namespace
