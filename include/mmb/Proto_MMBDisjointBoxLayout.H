#pragma once
#ifndef _PROTO_MMB_DISJOINT_BOX_LAYOUT_
#define _PROTO_MMB_DISJOINT_BOX_LAYOUT_
#include "Proto_DisjointBoxLayout.H"

namespace Proto
{

    class MMBDisjointBoxLayout
    {
    public:

    MMBDisjointBoxLayout();


    private:

    std::vector<DisjointBoxLayout> m_layouts;

    }; // end class MMBDisjointBoxLayout
}

#include "implem/Proto_MMBDisjointBoxLayoutImplem.H"
#endif // end include guard
