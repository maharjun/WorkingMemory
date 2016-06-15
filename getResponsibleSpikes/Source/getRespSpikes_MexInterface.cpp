#include <mex.h>
#include <matrix.h>
#undef printf

#include <MexMemoryInterfacing/Headers/GenericMexIO.hpp>
#include <MexMemoryInterfacing/Headers/InterruptHandling.hpp>

#include "../Headers/getRespSpikes.hpp"

#ifdef _MSC_VER
#  define SPRINTF_FUNC sprintf_s
#  define STRTOK_FUNC strtok_s
#  define STRCMPI_FUNC _strcmpi
#elif defined __GNUC__
#  if (__GNUC__ > 5) || (__GNUC__ == 5)
#    define SPRINTF_FUNC std::snprintf
#    define STRTOK_FUNC strtok_r
#    define STRCMPI_FUNC strcasecmp
#  endif
#endif

using namespace std;

void mexFunctionCpp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// This function is the function that does the actual work of the MEX
	// Function including performing IO from and to mxArrays, and calling the
	// Simulation function.

	// Validate Number of Arguments
	if (nrhs != 1) {
		WriteException(ExOps::EXCEPTION_INVALID_INPUT, "Incorrect number of input arguments (required 1, given %d)", nrhs);
	}
	if (nlhs != 1) {
		WriteException(ExOps::EXCEPTION_INVALID_INPUT, "Incorrect number of output arguments (required 1, given %d)", nlhs);
	}

	// Initialize Simulation Objects and do simulation
	getRespSpikes::SimulationClass SimulationObj;
	SimulationObj.initialize(prhs[0]);

	// Running Simulation
	chrono::system_clock::time_point TStart = chrono::system_clock::now();
	EnableInterruptHandling();
	SimulationObj.ResponsibleSynCalc();
	DisableInterruptHandling();
	chrono::system_clock::time_point TEnd = chrono::system_clock::now();

	WriteOutput("The Time taken = %d milliseconds\n", chrono::duration_cast<chrono::milliseconds>(TEnd - TStart).count());

	plhs[0] = SimulationObj.getOutput();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// This function is the C wrapper to mexFunctionCpp that performs
	// Memory account opening and exception handling. This function is only
	// called when compiled into a MEX file (i.e. MEX_LIB is defined)

	// Open Memory Usage Account
	size_t MemAccountKey = MemCounter::OpenMemAccount(size_t(4) << 29);

	try {
		mexFunctionCpp(nlhs, plhs, nrhs, prhs);
	}
	catch (ExOps::ExCodes A) {
		if (A == ExOps::EXCEPTION_MEM_FULL) {
			char OutputString[256];
			SPRINTF_FUNC(OutputString, 256, "EXCEPTION: Mem Limit of %lld MB Exceeded\n", (MemCounter::MemUsageLimit) >> 20);
			mexErrMsgIdAndTxt("CppSimException:MemOverFlow", OutputString);
		}
		else if (A == ExOps::EXCEPTION_INVALID_INPUT) {
			char OutputString[256];
			SPRINTF_FUNC(OutputString, 256, "EXCEPTION: Invalid Input\n");
			mexErrMsgIdAndTxt("CppSimException:InvalidInput", OutputString);
		}
		else if (A == ExOps::EXCEPTION_CONST_MOD || A == ExOps::EXCEPTION_EXTMEM_MOD) {
			char OutputString[256];
			SPRINTF_FUNC(OutputString, 256, "EXCEPTION: Invalid Modification of %s Memory\n", ((A == ExOps::EXCEPTION_CONST_MOD)?"const":"external") );
			mexErrMsgIdAndTxt("CppSimException:InvalidInput", OutputString);
		}
	}

	// Close Memory Usage Account
	MemCounter::CloseMemAccount(MemAccountKey);
}