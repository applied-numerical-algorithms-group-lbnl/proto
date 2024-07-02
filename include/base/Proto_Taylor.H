#ifndef _PROTO_TAYLOR_H_
#define _PROTO_TAYLOR_H_

#include <cassert>
#include <cmath>

namespace Proto {
    
    int factorial(int value)
    {
        if (value == 0) { return 1; }
        else { return value * factorial(value - 1); }
    }

    template<typename T>
    T taylor(int index) 
    {
        return 1.0/((T)factorial(index));
    }

} //end Proto namespace

#endif //end include guard
