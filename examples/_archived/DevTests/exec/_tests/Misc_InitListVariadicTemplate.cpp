#include <iostream>
#include <array>

void foo(std::initializer_list<std::array<int, 3>> a_args)
{
    for (auto iter : a_args)
    {
        for (int ii = 0; ii < 3; ii++)
        {
            std::cout << iter[ii] << ", ";
        }
        std::cout << std::endl;

    }
}
int main(int argc, char** argv)
{
    foo({{1,2,3}, {4,5,6}});
}
