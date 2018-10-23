#include <cstring>

extern "C" {
inline int omp_get_initial_device(){ return 0;}
inline void* omp_target_alloc(size_t size, int device){ return malloc(size);}
inline int omp_target_memcpy(void * dst, void * src, size_t length, size_t dst_offset,
                  size_t src_offset, int dst_device_num, int src_device_num)
{  memcpy((char*)dst+dst_offset, (char*)src+src_offset, length);
  return 0;
}
inline void omp_target_free(void * device_ptr, int device_num) {free(device_ptr);}
}
