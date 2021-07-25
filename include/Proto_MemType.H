<<<<<<< HEAD
#ifndef _Proto_MEMTYPE_H_
#define _Proto_MEMTYPE_H_
namespace Proto
{
enum MemType {HOST,DEVICE};
}
=======
#ifndef _PROTO_MEMTYPE_H_
#define _PROTO_MEMTYPE_H_

namespace Proto
{
  // enum class MemType {HOST,DEVICE};
  enum MemType {INVALID=0,HOST=1, DEVICE=2,BOTH=3};
#ifndef MEMTYPE_DEFAULT
#ifdef PROTO_CUDA
#define MEMTYPE_DEFAULT MemType::DEVICE
#else
#define MEMTYPE_DEFAULT MemType::HOST
#endif
#endif
}

>>>>>>> master
#endif
