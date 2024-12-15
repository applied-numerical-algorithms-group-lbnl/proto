#pragma once
#ifndef _PROTO_PATCH_LAYOUT_
#define _PROTO_PATCH_LAYOUT_

#include "Proto.H"

namespace Proto {

    class PatchLayout {
        public:
        PatchLayout(Domain patchDomain);
        PatchLayout(
            Domain patchDomain,
            std::vector<PointInfo> patches);
        void define(
            Domain patchDomain,
            std::vector<PointInfo> patches);

        inline void loadBalance(int startProc, int endProc);
        inline bool compatible(const PatchLayout& patchLayout) const;
        inline int numProcs() const;
        inline int numLocalPatches() const;
        inline int numGlobalPatches() const;
        inline int numPatchesOnProc(int proc) const;

        template<typename... args>
        inline void loadAssign(args... vals);

        private: 
        Domain m_patchDomain;
        std::vector<PatchID> m_patches;
        std::vector<std::pair<unsigned int, unsigned int>> m_procAssignment;

    };
}

#endif