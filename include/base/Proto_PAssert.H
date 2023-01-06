#pragma once
#ifndef __PROTO_ASSERT__
#define __PROTO_ASSERT__
#include <cassert>
#include <cstdio>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>

namespace ProtoUtils
{
    inline std::vector<std::string> split_string (const std::string &s, char delim) {
        std::vector<std::string> result;
        std::stringstream ss (s);
        std::string item;

        while (getline (ss, item, delim)) {
            result.push_back (item);
        }

        return result;
    }
}
// Custom Assert macro with printf-style reporting
#ifdef __CUDA_ARCH__

#define PROTO_ASSERT(stmt,args...)
#define PR_DEBUG_MSG(N, stmt, args...)
#define PR_assert(stmt) 
#define PR_error(stmt) 
#define PR_warning(stmt) 

#elif NDEBUG

#define PROTO_ASSERT(stmt,args...)
#define PR_DEBUG_MSG(N, stmt, args...)
#define PR_assert(stmt)
#define PR_error(stmt)
#define PR_warning(stmt)

#else


#ifndef PROTO_HIP

#define PROTO_ASSERT(stmt,args...) \
  assert( (stmt) ||                                         \
     (fprintf(stderr,"%s:%d: error: ",__FILE__,__LINE__) && \
     fprintf(stderr,"Assertion `" #stmt "` failed.\n\t") && \
     fprintf(stderr,args) && \
     fprintf(stderr,"\n"))); \
  if (!(stmt)){std::abort();}

#define PR_DEBUG_MSG(N, stmt, args...) \
    if (N <= PR_VERBOSE) \
    { \
        auto words = ProtoUtils::split_string(__FILE__, '/'); \
        string filename = words[words.size()-1]; \
        char prefix[100]; \
        char message[500]; \
        sprintf(prefix, "DEBUG | %s:%d: ", filename.c_str(), __LINE__); \
        sprintf(message, stmt, args); \
        fprintf(stdout,"%s %s\n", prefix, message); \
    } \

#define PR_assert(stmt) assert( (stmt) );


#define PR_error(stmt) \
  fprintf(stderr,"error thrown =  `" #stmt "`!!!\n\t"); \
  std::abort();

#define PR_warning(stmt) \
  fprintf(stderr,"warning thrown =  `" #stmt "` !!!\n\t");

#else // PROTO_HIP
#define PROTO_ASSERT(stmt,args...)
#define PR_assert(stmt) 
#define PR_error(stmt) 
#define PR_warning(stmt) 

#endif
#endif
#endif // end include guard
