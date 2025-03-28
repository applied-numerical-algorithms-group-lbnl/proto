#pragma once
#include <iostream>

namespace Classes {
    class A {
        private:
        int value;
        
        public:

        inline A();
        inline void print();

    };
    #include "implem/ClassAImplem.H"
}
