#pragma once
#ifndef _PROTO_BOX_LAYOUT_
#define _PROTO_BOX_LAYOUT_

#include "Proto.H"

namespace Proto {

    class BoxLayout {
        
        public:

        inline BoxLayout(Domain domain);
        inline BoxLayout(Domain domain, std::vector<PatchID> patches);
        inline void define(Domain domain);
        inline void define(Domain domain, std::vector<PatchID> patches);
        inline void define(std::shared_ptr<PatchLayout> patchLayout, Point boxSize);
        inline void define(std::shared_ptr<PatchLayout> patchLayout, std::vector<Point> boxSizes);

        private:
        Domain m_domain;
        std::shared_ptr<PatchLayout> m_patchLayout;
    };
}
#endif