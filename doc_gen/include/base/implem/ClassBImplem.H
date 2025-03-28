#pragma once
#include "../ClassB.H"
using namespace Classes;

B::B() { values = {1,2,3}; }
void B::print()
{
    for (auto vi : values)
    {
        std::cout << vi << ", ";
    }
    std::cout << std::endl;
}