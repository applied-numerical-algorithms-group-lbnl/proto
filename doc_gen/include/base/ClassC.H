#pragma once
#include "ClassA.H"
#include "ClassB.H"

namespace Classes {
    
    typedef A typedefA;
    typedef B typedefB;
    
    class C {
        private: 

        A a;
        B b;

        public:

        void print();
    };
}