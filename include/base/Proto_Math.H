#ifndef _PROTO_MATH_H_
#define _PROTO_MATH_H_

namespace Proto {
    
    /// Template Based Integer Exponentiation
    template <unsigned int P>
    inline int ipow(int M){return M*ipow<P-1>(M);}
    
    template <>
    inline int ipow<0>(int M){return 1;}

    // Dynamic version
    inline int ipow(int a_base, unsigned int a_exp)
    {
        int result = 1;
        for(;;)
        {
            if (a_exp & 1){result *= a_base;}
            a_exp >>= 1;
            if (!a_exp){break;}
            a_base *= a_base;
        }
        return result;
    }

    inline int factorial(int value)
    {
        if (value == 0) { return 1; }
        else { return value * factorial(value - 1); }
    }

    template<typename T>
    inline T taylor(int index) 
    {
        return 1.0/((T)factorial(index));
    }

} //end Proto namespace

#endif //end include guard
