#include <iostream>

enum Value
{
    A,
    B,
    C
};

template<Value V1, Value V2>
void foo()
{
#if V1 == 1
#if V2 == 1
    std::cout << "Values are both B" << std::endl;
#endif
#endif
};

int main(int argc, char** argv)
{
    foo<B,C>();
    foo<A,B>();
    foo<B,B>();
};
