#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    for (auto iter = Side::begin(); iter != Side::end(); ++iter)
    {
        std::cout << *iter << std::endl;
    }
}
