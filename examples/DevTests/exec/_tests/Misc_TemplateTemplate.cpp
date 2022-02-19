#include <iostream>


template<unsigned int C>
class A
{
    public:
    inline static constexpr unsigned int N() {return C;}
    A() {m_value = C;}
    void print() {std::cout << m_value << std::endl; }

    private:
    int m_value;
};

class A7 : public A<7>
{
};

template<class T>
class B
{
    public:

    typedef A<T::N()> V;

    B(){}
    void print(){ std::cout << V::N() << std::endl; }
};

int main(int argc, char** argv)
{
    B<A7> b;
    b.print();
}
