#include "Proto.H"

using namespace Proto;

template<typename T>
struct getMemType_
{
    static constexpr MemType type()
    {
        return Proto::MemType::BOTH;
    }
    static MemType type_eval()
    {
        return Proto::MemType::BOTH;
    }
};

template<typename T, unsigned int C, MemType MEM, unsigned char D, unsigned char E>
struct getMemType_<BoxData<T, C, MEM, D, E>>
{
    static constexpr MemType type() 
    {
        return MEM;
    }
    static MemType type_eval()
    {
        return MEM;
    }
};


int main(int argc, char** argv)
{
    BoxData<double, 1, HOST> hostData;
    BoxData<double, 1, DEVICE> deviceData;

    std::cout << "General value: " << getMemType_<double>::type_eval() << std::endl;
    std::cout << "Host BoxData value: " << getMemType_<decltype(hostData)>::type_eval() << std::endl;
    std::cout << "Device BoxData value: " << getMemType_<decltype(deviceData)>::type_eval() << std::endl;
    std::cout << "Host BoxData value: " << getMemType<decltype(hostData)>::type_eval() << std::endl;
    std::cout << "Device BoxData value: " << getMemType<decltype(deviceData)>::type_eval() << std::endl;
    
}




