#pragma once
#include "Proto_BoxLayout.H"

namespace Proto {

// forward declarations
class LevelCopier;

template <class T, unsigned int C, unsigned char D, unsigned char E>
class LevelData
{
public:
    inline void define(const BoxLayout& a_layout,
                       Point a_ghost = Point::Zeros())
    {
        //set m_rank via MPI / UPC++
        auto iter = a_layout.dataIterator();
        for (iter.begin(); !iter.done(); ++iter)
        {
            BoxData<T,C,D,E> patch(a_layout[iter].grow(a_ghost));
            m_patches[*iter] = &patch;
        }
    }

    inline void exchange();
    inline void copyTo(LevelData<T,C,D,E>& a_dest) const;

    inline BoxData<T,C,D,E>& operator[](const DataIterator& a_iter)
    {
        return &m_patches[*a_iter];
    }

     
private:
    int m_rank;                                 ///< This parallel rank
    BoxLayout m_layout;                         ///< Proto::BoxLayout equivalent to Chombo's DisjointBoxLayout. 
    std::map<Point,BoxData<T,C,D,E>*> m_patches; ///< BoxData patches local to this proc
    LevelCopier<T,C,D,E> m_copier;              ///< To facilitate copyTo / exchange
}

template<class T, unsigned int C, unsigned char D, unsigned char E>
class LevelCopier
{

}

} //end Proto namespace
