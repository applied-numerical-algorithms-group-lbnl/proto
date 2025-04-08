#pragma once
#include "../ClassB.H"

namespace Classes {

    inline B::B() { values = {1,2,3}; }
    inline void B::print()
    {
        for (auto vi : values)
        {
            std::cout << vi << ", ";
        }
        std::cout << std::endl;
    }
}