#include "InputParser.H"

int main(int argc, char** argv)
{
    InputArgs args;
    args.parse(argc, argv);
    args.print();
};
