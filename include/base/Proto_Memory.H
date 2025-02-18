#pragma once
#ifndef _PROTO_MEMORY_
#define _PROTO_MEMORY_

#include <cstring>
#include <limits>
#include "Proto_MemType.H"

#include <sys/resource.h>

#define PR_UNALLOCATED_FLOAT 7.77e77

//#define   RUSAGE_SELF     0
//#define   RUSAGE_CHILDREN     -1

namespace Proto
{

// The purpose of these functions is to improve the existing macro-driven
// API for memory access in Proto. Existing macros have not been removed,
// but will be phased out in favor of templated functions as development
// progresses.

// TODO: TECHNICALLY the MemType template parameters might be redundant. 
// Both HIP and CUDA have an API for checking of a pointer is a device 
// pointer or not, but based on my limited testing, I'm not sure how
// much faith I have in said tools. 

/// Copy Memory
/**
    Executes a low level memory copy between buffers of the specified MemTypes.
    If both SRC_MEM and DST_MEM are HOST, this is equivalent to <code>std::memcpy</code>.
    When at least one MemType is DEVICE, this function will call 
    <code>cudaMemcpy</code> or equivalent functions as appropriate.

    \tparam SRC_MEM     MemType of the source buffer
    \tparam DST_MEM     MemType of the destination buffer
    \param a_src        Source buffer
    \param a_dst        Destination buffer
    \param a_nbytes     Number of bytes to copy
*/
template<MemType SRC_MEM=MEMTYPE_DEFAULT, MemType DST_MEM=MEMTYPE_DEFAULT>
inline void proto_memcpy(const void* a_src, void* a_dst, size_t a_nbytes);

/// Allocate Memory
/**
    Returns a pointer to a region of allocated memory. If <code>MEM==HOST</code>
    this function is equivalent to <code>std::malloc</code>.
    If <code>MEM==DEVICE</code>, this function will call <code>cudaMalloc</code>
    or other equivalent functions as appropriate.

    \tparam MEM     MemType of the desired buffer
    \param a_nbytes Number of bytes to allocate
*/
template<MemType MEM=MEMTYPE_DEFAULT>
inline void* proto_malloc(size_t a_nbytes);

/// Free Memory
/**
    Frees a pointer previously allocated using <code>proto_malloc</code>.
    If <code>MEM==HOST</code> this function is equivalent to <code>std::free</code>.
    If <code>MEM==DEVICE</code>, this function will call <code>cudaFree</code>
    or other equivalent functions as appropriate.

    \tparam MEM     MemType of the buffer.
    \param a_buffer A buffer.
*/
template<MemType MEM=MEMTYPE_DEFAULT>
inline void proto_free(void* a_buffer);

/// Query Pointer MemType
/**
    This is an EXPERIMENTAL tool for checking the MemType of an arbitrary pointer.
    Initial testing suggests that this function can distinguish between HOST and
    DEVICE for cuda pointers. Untested with HIP. Use with healthy skepticism.

    \param a_ptr    A pointer.
*/
inline MemType pointerMemType(const void* a_ptr);

template<typename T>
inline constexpr T proto_unallocated()
{
    if constexpr (std::is_floating_point_v<T>)
    {
        T val = (T)PR_UNALLOCATED_FLOAT;
        return val;
    } else {
        return std::numeric_limits<T>::max()-1;
    }
}

#include "implem/Proto_MemoryImplem.H"
}
#endif
