#ifndef __DEBUGHOOKS__
#define __DEBUGHOOKS__

#include "WriteBoxData.H"

void printData1(BoxData<double, 1>* dataPtr);
void printDataDim(BoxData<double, DIM>* dataPtr);
void writeData1(BoxData<double, 1>* dataPtr);
void maxMinData1(BoxData<double, 1>* dataPtr);
void maxMinDataDim(BoxData<double, DIM>* dataPtr);
void  viewData1(BoxData<double, 1>* dataPtr);
void writeDataDim(BoxData<double, DIM>* dataPtr);
void  viewDataDim(BoxData<double, DIM>* dataPtr);
void ProtoFunctionToGetStuffIntoTheSymbolTable();

             
             

#include "DebugHooksImplem.H"
#endif
