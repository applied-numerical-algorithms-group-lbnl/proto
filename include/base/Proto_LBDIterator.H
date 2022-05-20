#pragma once
#ifndef _PROTO_LBD_ITERATOR_
#define _PROTO_LBD_ITERATOR_

namespace Proto {

    template<typename T, unsigned int C, MemType MEM>
    class LBDIterator
    {
        public:
            LBDIterator(LevelBoxData<T, C, MEM>& a_data);
            BoxData<T, C, MEM>& operator*();
        private:
            LevelBoxData<T, C, MEM>& m_data;
            DataIndex m_index;
    };

} // end namespace Proto

#endif // end include guard
