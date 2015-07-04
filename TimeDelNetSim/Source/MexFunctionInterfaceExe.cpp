#include <matrix.h>
#include <mat.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <chrono>
#include <type_traits>
#include <iostream>
#include <Windows.h>
#include "..\Headers\NeuronSim.hpp"
#include "..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "MexFunctionInterface.cpp"

using namespace std;

typedef mxArray* mxArrayPtr;

int main(){
	// NOTE THAT THERE IS NO DATA VALIDATION AS THIS IS EXPECTED TO HAVE 
	// BEEN DONE IN THE MATLAB SIDE OF THE INTERFACE TO THIS MEX FUNCTION
	mxArrayPtr  Input      = nullptr, 
				InitState  = nullptr,
				OutputVars = nullptr,
				StateVars  = nullptr,
				FinalState = nullptr;

	MATFile* InputFilePtr = matOpen("Data/InputData.mat", "r");
	MATFile* OutputFilePtr = nullptr;
	char OutFileName[256];
	char OutputFilePath[256] = "";
	Input = matGetVariable(InputFilePtr, "InputStruct");
	
	mxGetString_730(mxGetField(Input, 0, "OutputFile"), OutFileName, 256);
	strcat_s(OutputFilePath, 256, "Data/");
	strcat_s(OutputFilePath, 256, OutFileName);

	matClose(InputFilePtr);
	if (Input == nullptr){
		cout << "TimeDelNetSimEXE:InvInpVarName", "The variable name in the mex file InputData must be InputStruct";
	}
	
	mxArrayPtr lhs[4] = { nullptr, nullptr, nullptr, nullptr }, 
			   rhs[1] = { Input };

	HANDLE XHandle = GetCurrentProcess();
	ULONG_PTR ProcessAffMask, SystemAffMask;
	
	GetProcessAffinityMask(XHandle, &ProcessAffMask, &SystemAffMask);
	cout << std::hex << ProcessAffMask << "  " << std::hex << SystemAffMask << endl;

	system("pause");

	OutputFilePtr = matOpen(OutputFilePath, "r");
	while (OutputFilePtr){
		char UserConfirmResp;
		std::cout << "File Exists. Sure about rewrite? : ";
		std::cin >> UserConfirmResp;
		if ((UserConfirmResp | 32) == 'y'){
			matClose(OutputFilePtr);
			OutputFilePtr = nullptr;
		}
		else if ((UserConfirmResp | 32) == 'n'){
			matClose(OutputFilePtr);
			OutputFilePtr = nullptr;
			cout << "KTHXBYE" << endl;
			system("pause");
			mxDestroyArray(Input);
			return 0;
		}
	}
	try{
		mexFunction(4, lhs, 1, rhs);
	}
	catch (ExOps::ExCodes A){
		if (A == ExOps::EXCEPTION_MEM_FULL){
			printf("Mem Limit of %lld MB Exceeded\n", (MemCounter::MemUsageLimit) >> 20);
			system("pause");
		}
		mxDestroyArray(Input);
		return 0;
	}
	
	OutputVars = lhs[0];
	StateVars = lhs[1];
	FinalState = lhs[2];
	InitState = lhs[3];
	
	OutputFilePtr = matOpen(OutputFilePath, "wz");
	matPutVariable(OutputFilePtr, "InitState", InitState);
	matPutVariable(OutputFilePtr, "OutputVars", OutputVars);
	matPutVariable(OutputFilePtr, "StateVars", StateVars);
	matPutVariable(OutputFilePtr, "FinalState", FinalState);
	matClose(OutputFilePtr);
	OutputFilePtr = nullptr;
	
	mxDestroyArray(Input);
	mxDestroyArray(InitState);
	mxDestroyArray(OutputVars);
	mxDestroyArray(StateVars);
	mxDestroyArray(FinalState);
	system("pause");
	return 0;
}