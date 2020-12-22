#ifndef _PROTO_MEMTYPE_H_
#define _PROTO_MEMTYPE_H_

namespace Proto
{
    enum class MemType {HOST,DEVICE};
#ifndef MEMTYPE_DEFAULT
#ifdef PROTO_CUDA
#define MEMTYPE_DEFAULT MemType::DEVICE
#else
#define MEMTYPE_DEFAULT MemType::HOST
#endif
#endif
}

#endif
