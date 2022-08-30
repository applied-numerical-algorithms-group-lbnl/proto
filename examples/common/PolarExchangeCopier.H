#pragma once
#ifndef _POLAR_MIRROR_COPIER_H_
#define _POLAR_MIRROR_COPIER_H_

#include "Proto.H"

using namespace Proto;

/// Polar Mirror Copier
/** A custom Copier which executes a LevelData exchange operation at the pole of a grid
 * represented in spherical coordinates. The behavior is identical to the 
 * default ExchangeCopier except at boundaries in the elevation coordinate direction
 * in which a boundary cells are filled with data shifted by 180* in the azimuthal direction
*/
template<typename T, unsigned int C, MemType MEM, Centering CTR>
class PolarExchangeCopier
: public LevelExchangeCopier<T, C, MEM, CTR>
{
    public:

    /// Constructor
    /**
     * \param a_azimuthal   The azimuthal coordinate (0 - 360)
     * \param a_elevation   The elevation coordinate (0 - 180)
    */
    inline PolarExchangeCopier(unsigned char a_azimuthal, unsigned char a_elevation)
    {
        m_azimuthal = a_azimuthal;
        m_elevation = a_elevation;
    };

    /// Build Copier Motion Plan
    inline void buildMotionPlans(LevelCopierOp<T, C, MEM, MEM, CTR>& a_op);

    private:
    unsigned char m_azimuthal;
    unsigned char m_elevation;
};


template<typename T, unsigned int C, MemType MEM, Centering CTR>
void PolarExchangeCopier<T, C, MEM, CTR>::buildMotionPlans(LevelCopierOp<T, C, MEM, MEM, CTR>& a_op)
{
    auto& data = *a_op.m_src;
    const auto& layout = data.layout();
    Box domainBox = layout.domain().box();
    Box patchDomainBox = layout.patchDomain().box();
    PROTO_ASSERT(domainBox.size(m_azimuthal) % 2 == 0,
            "PolarExchangeCopier::buildMotionPlans | Error: domain box must have an even number of cells \
            in the azimuthal direction.");
    PROTO_ASSERT(patchDomainBox.size(m_azimuthal) % 2 == 0,
            "PolarExchangeCopier::buildMotionPlans | Error: layout must have an even number of patches \
            in the azimuthal direction.");

    LevelExchangeCopier<T, C, MEM, CTR>::buildMotionPlans(a_op);

    Point ghost = data.ghost();
    Box ghostDomainBox = domainBox.grow(ghost);
    int ghostSize = ghost[m_elevation];
    Point e_theta = Point::Basis(m_azimuthal);
    Point e_phi = Point::Basis(m_elevation);
    Box domainNorth = ghostDomainBox.edge(e_phi, ghostSize);
    Box domainSouth = ghostDomainBox.edge(-e_phi, ghostSize);
    Box bdryPatchesNorth = patchDomainBox.grow(1).edge(e_phi,1);
    Box bdryPatchesSouth = patchDomainBox.grow(1).edge(-e_phi,1);
    Point mirrorShift = e_theta*domainBox.size(m_azimuthal)/2; 
    Point mirrorPatchShift = mirrorShift / layout.boxSize()[m_azimuthal]; 

    CoordPermutation mirror({{m_elevation, m_elevation, -1}});

    auto& localPlan = this->m_localMotionPlan;
    auto liter = localPlan.begin();
    while(liter != localPlan.end())
    {
        const Box& toRegion = (*liter).toRegion;
        if (domainNorth.contains(toRegion) || domainSouth.contains(toRegion))
        {
            liter = localPlan.erase(liter);
        } else {
            ++liter;
        }
    }
    auto& fromPlan = this->m_fromMotionPlan;
    auto fiter = fromPlan.begin();
    while(fiter != fromPlan.end())
    {
        const Box& toRegion = (*fiter).toRegion;
        if (domainNorth.contains(toRegion) || domainSouth.contains(toRegion))
        {
            fiter = fromPlan.erase(fiter);
        } else {
            ++fiter;
        }
    }
    auto& toPlan = this->m_toMotionPlan;
    auto titer = toPlan.begin();
    while(titer != toPlan.end())
    {
        const Box& toRegion = (*titer).toRegion;
        if (domainNorth.contains(toRegion) || domainSouth.contains(toRegion))
        {
            titer = toPlan.erase(titer);
        } else {
            ++titer;
        }
    }

    Box dirs = Box::Kernel(1);
    for (auto iter : layout)
    {
        Point toPatch = layout.point(iter);
        for (auto dir : dirs)
        {
            Point n = toPatch+dir;
            if (bdryPatchesNorth.contains(n) || bdryPatchesSouth.contains(n))
            {
                Point norm = dir*e_phi;
                Point fromPatch = n - norm + mirrorPatchShift;
                Point reverseDir = -(dir-norm)+norm;
                auto fromIndex = layout.find(fromPatch);
                if (fromIndex != layout.end())
                {
                    Box toHereBox = layout[iter].adjacent(dir*ghost);
                    Box fromThereBox = layout[fromIndex].edge(reverseDir*ghost);
                    LevelMotionItem toItem(fromIndex, iter, fromThereBox, toHereBox);
                    toItem.rotation = mirror;
                    auto fromProc = layout.procID(fromIndex);
                    if (fromProc == procID())
                    {
                        this->m_localMotionPlan.push_back(toItem);
                    } else {
                        toItem.procID = fromProc;
                        this->m_toMotionPlan.push_back(toItem);
                        
                        Box toThereBox = layout[fromIndex].adjacent(reverseDir*ghost);
                        Box fromHereBox = layout[iter].edge(dir*ghost);
                        LevelMotionItem fromItem(
                                iter, fromIndex, fromHereBox, toThereBox);
                        fromItem.rotation = mirror;
                        fromItem.procID = fromProc;
                        this->m_fromMotionPlan.push_back(fromItem);
                    }
                }
            }
        }
    }
    this->sort();
}
#endif //end include guard
