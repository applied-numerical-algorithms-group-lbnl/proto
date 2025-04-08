#pragma once
#include <iostream>

namespace Classes {
    
    enum AEnum {
        E0,
        E1,
        E2
    };
    
    class A {
        private:
        int value;
        
        public:

        A();
        void print() const;

    };

    void fooA(const A& a);
}
