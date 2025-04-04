#pragma once
#include "../ClassA.H"
namespace Classes {
    inline A::A() { value = 7; }
    inline void A::print() const { std::cout << value << std::endl; }
    void fooA(const A& a)
    {
        a.print();
    }
}
