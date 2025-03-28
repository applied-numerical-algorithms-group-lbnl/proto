#pragma once
#include <iostream>
#include <vector>

namespace Classes {
    class B {
        private:
        std::vector<int> values;
        
        public:

        inline B();
        inline void print();

    };
    #include "implem/ClassBImplem.H"
}