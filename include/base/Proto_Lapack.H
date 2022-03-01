#pragma once
#ifndef _PROTO_LAPACK_
#define _PROTO_LAPACK_

extern "C" {
    int dgemv_(char*, int*, int*, double*, double*,
            int*, double*, int*, double*, double*, int*);
    int dgemm_(char*, char*, int*, int*, int*, double*,
            double*, int*, double*, int*, double*, double*, int*);
    int dgetri_(int*, double*, int*, int*, double*, int*, int*);
    int dgetrf_(int*, int*, double*, int*, int*, int*);
}

#endif // end include guard
