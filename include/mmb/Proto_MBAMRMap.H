#pragma once
#ifndef _PROTO_MB_AMR_MAP_
#define _PROTO_MB_AMR_MAP_

#include "Proto_MBAMRData.H"

namespace Proto
{

template<MemType MEM = MEMTYPE_DEFAULT>
class MBAMRMap
{

};
#include "implem/Proto_MBAMRMapImplem.H"
} //end Proto namespace
