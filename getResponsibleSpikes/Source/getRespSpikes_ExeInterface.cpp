#include <matrix.h>
#include <mat.h>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <chrono>
#include <type_traits>
#include <iostream>

#include <MexMemoryInterfacing/Headers/MexMem.hpp>
#include <MexMemoryInterfacing/Headers/GenericMexIO.hpp>

#ifdef _MSC_VER
#  define STRCAT_SAFE(a,b,c) strcat_s((a),(b),(c))
#elif defined __GNUC__
#  if (__GNUC__ > 5) || (__GNUC__ == 5)
#    define STRCAT_SAFE(a,b,c) strncat((a),(c),(b))
#  endif
#endif

extern void mexFunctionCpp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

using namespace std;

int main(int argc, char *args[]){

	// Open Memory Usage Account
	size_t MemAccountKey =  MemCounter::OpenMemAccount(size_t(4) << 29);

	// Start Data Input from MAT File
	mxArrayPtr Input = nullptr, RespSpikesStruct = nullptr;

	char InputFilePath[256] = "";
	char OutFileName[256] = "";
	char OutputFilePath[256] = "";

	// Process Input Aruments and get filenames
	if (argc == 1) {
		STRCAT_SAFE(InputFilePath, 256, "InputData.mat");
		STRCAT_SAFE(OutputFilePath, 256, "OutputData.mat");
	}
	else if(argc == 2) {
		STRCAT_SAFE(InputFilePath, 256, args[1]);
		STRCAT_SAFE(OutputFilePath, 256, "OutputData.mat");
	}
	else if (argc == 3) {
		STRCAT_SAFE(InputFilePath, 256, args[1]);
		STRCAT_SAFE(OutputFilePath, 256, args[2]);
	}
	else {
		WriteOutput("Require at max two arguments\n");
		std::cin.ignore(1024, '\n');
		std::cout << "Press enter to terminate simulation...";
		std::cin.get();
	}

	MATFile* InputFilePtr = matOpen(InputFilePath, "r");
	MATFile* OutputFilePtr = nullptr;

	if (InputFilePtr == nullptr) {
		WriteException(ExOps::EXCEPTION_INVALID_INPUT, "The File '%s' does not exist.\n", InputFilePath);
	}

	Input = matGetVariable(InputFilePtr, "InputStruct");
	if (Input == nullptr){
		WriteException(ExOps::EXCEPTION_INVALID_INPUT, "The variable name in the Input Data file must be 'InputStruct'");
	}

	mxArrayPtr lhs[1] = { nullptr };
	const mxArray* rhs[1] = { Input };

	// Confirm Output File Rewrite
	OutputFilePtr = matOpen(OutputFilePath, "r");
	while (OutputFilePtr){
		char UserConfirmResp;
		std::cout << "The following file already exists - \n" << std::endl << "    " << OutputFilePath << std::endl << "\nSure about rewrite? : " << std::flush;
		std::cin >> UserConfirmResp;
		if ((UserConfirmResp | 32) == 'y'){
			matClose(OutputFilePtr);
			OutputFilePtr = nullptr;
		}
		else if ((UserConfirmResp | 32) == 'n'){
			matClose(OutputFilePtr);
			OutputFilePtr = nullptr;
			cout << "KTHXBYE" << endl;

			std::cin.ignore(1024, '\n');
			std::cout << "Press enter to terminate simulation...";
			std::cin.get();

			mxDestroyArray(Input);
			return 0;
		}
	}

	// Execute Simulation via the Mex Interface
	try{
		mexFunctionCpp(1, lhs, 1, rhs);
	}
	catch (ExOps::ExCodes A){
		if (A == ExOps::EXCEPTION_MEM_FULL){
			WriteOutput("Mem Limit of %lld MB Exceeded\n", (MemCounter::MemUsageLimit) >> 20);
		}
		else if (A == ExOps::EXCEPTION_INVALID_INPUT) {
			WriteOutput("EXCEPTION: Invalid Input\n");
		}
		else if (A == ExOps::EXCEPTION_CONST_MOD || A == ExOps::EXCEPTION_EXTMEM_MOD) {
			WriteOutput("Invalid Modification of %s Memory\n", ((A == ExOps::EXCEPTION_CONST_MOD) ? "const" : "external"));
		}

		mxDestroyArray(Input);
		for (int j = 0; j < 4; ++j) {
			if (lhs[j] != nullptr)
				mxDestroyArray(lhs[j]);
		}

		std::cin.ignore(1024, '\n');
		std::cout << "Press enter to terminate simulation...";
		std::cin.get();

		return 0;
	}

	RespSpikesStruct = lhs[0];

	// Write Output into Output MAT File
	OutputFilePtr = matOpen(OutputFilePath, "wz");
	matPutVariable(OutputFilePtr, "RespSpikesStruct", RespSpikesStruct);
	matClose(OutputFilePtr);
	OutputFilePtr = nullptr;

	mxDestroyArray(Input);
	mxDestroyArray(RespSpikesStruct);

	// Close Memory Usage Account
	MemCounter::CloseMemAccount(MemAccountKey);

	std::cin.ignore(1024, '\n');
	std::cout << "Press enter to end simulation...";
	std::cin.get();
	return 0;
}