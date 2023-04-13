#include "InputParser.H"

int main(int argc, char** argv)
{

    double D = 1.7e-3;
    int I = 7;
    bool B = true;
    std::string S = "message";

    InputArgs args;
    args.add("D", D);
    args.add("I", I);
    args.add("B", B);
    args.add("S", S);
    
    args.parse(argc, argv);
    args.print();

    pout() << "D: " << D << std::endl;
    pout() << "I: " << I << std::endl;
    pout() << "B: " << B << std::endl;
    pout() << "S: " << S << std::endl;

};
