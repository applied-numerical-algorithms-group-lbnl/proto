#pragma once
#include "ClassA.H"
#include "ClassB.H"

namespace Classes {
    class C {
        private: 

        A a;
        B b;

        public:

        inline void print();
    };
    #include "implem/ClassCImplem.H"
}