#pragma once
#include <iostream>

namespace Classes {
    class A {
        private:
        int value;
        
        public:

        A();
        void print() const;

    };

    void fooA(const A& a);
}
