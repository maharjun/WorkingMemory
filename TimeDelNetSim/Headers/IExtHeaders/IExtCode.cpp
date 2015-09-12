#include "IExtCode.hpp"
#include <cstring>
#include <cctype>
#include <string>
#include <algorithm>

#include "..\Network.hpp"
#include "..\NeuronSim.hpp"

#include <matrix.h>
#include "..\..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\..\..\MexMemoryInterfacing\Headers\GenericMexIO.hpp"
#include "..\..\..\RandomNumGen\Headers\FiltRandomTBB.hpp"

////////////////////////////////////////////////////////
// Input and Initialization functions 
////////////////////////////////////////////////////////
size_t IExtInterface::getOutputControl(char * OutputControlString)
{
	MexVector<std::string> OutputControlOptions;
	StringSplit(OutputControlString, " ,-", OutputControlOptions);
	
	// Defining case insensitive comparison function
	auto CaseInsensitiveCharComp = [](char c1, char c2) -> bool { return std::tolower(c1) == std::tolower(c2);};
	auto iEqual = [&](const std::string &s1, const std::string &s2) -> bool {
		return std::equal(s1.begin(), s1.end(), s2.begin(), s2.end(), CaseInsensitiveCharComp);
	};

	// Defining return variable
	size_t OutputControlWord = 0;

	for (auto OutContOpt : OutputControlOptions) {
		// Split the current OutputControlOption by '.'
		MexVector<std::string> OutputControlOptionParts;
		StringSplit(OutContOpt.data(), ".", OutputControlOptionParts);

		bool AddorRemove = true; // TRUE for Add

		if (iEqual(OutputControlOptionParts[0], "FSF")) {
			OutputControlWord |=
				IExtInterface::OutOps::GEN_STATE_REQ |
				IExtInterface::OutOps::I_RAND_REQ;
		}

		if (iEqual(OutputControlOptionParts[0], "IExt") || iEqual(OutputControlOptionParts[0], "/IExt")
			&& OutputControlOptionParts.size() == 2) {
			// Ascertain Add or Remove
			if (OutputControlOptionParts[0][0] == '/') {
				AddorRemove = false;
			}

			// Check out different Output Option cases
			if (iEqual(OutputControlOptionParts[1], "GenState")) { 
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::GEN_STATE_REQ : 
					OutputControlWord & ~IExtInterface::OutOps::GEN_STATE_REQ;
			}
			if (iEqual(OutputControlOptionParts[1], "Irand")) {
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_RAND_REQ : 
					OutputControlWord & ~IExtInterface::OutOps::I_RAND_REQ;
			}
			if (iEqual(OutputControlOptionParts[1], "IextWORand")) {
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_EXT_WO_RAND_REQ : 
					OutputControlWord & ~IExtInterface::OutOps::I_EXT_WO_RAND_REQ;
			}
			if (iEqual(OutputControlOptionParts[1], "IextTotal")) {
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_EXT_TOTAL_REQ : 
					OutputControlWord & ~IExtInterface::OutOps::I_EXT_TOTAL_REQ;
			}
		}
	}

	return OutputControlWord;
}

void IExtInterface::takeInputVarsFromMatlabStruct(
	IExtInterface::InputVarsStruct & IExtInputVarsStruct,
	mxArray * IExtMatlabInputStruct, 
	InputArgs & SimulationInputArgs)
{
	// Giving Default Values to Optional Simulation Algorithm Parameters
	IExtInputVarsStruct.alpha = 0.5f;
	IExtInputVarsStruct.StdDev = 3.5f;

	// Taking input for Optional Simulation Algorithm Parameters
	getInputfromStruct<float>(IExtMatlabInputStruct, "Iext.alpha" , IExtInputVarsStruct.alpha );
	getInputfromStruct<float>(IExtMatlabInputStruct, "Iext.StdDev", IExtInputVarsStruct.StdDev);

	// Initializing OutputControl
	// Get OutputControlString and OutputControl Word
	mxArrayPtr genmxArrayPtr = getValidStructField(IExtMatlabInputStruct, "OutputControl", MexMemInputOps());
	if (genmxArrayPtr != NULL && !mxIsEmpty(genmxArrayPtr)) {
		char * OutputControlSequence = mxArrayToString(genmxArrayPtr);
		IExtInputVarsStruct.OutputControl = IExtInterface::getOutputControl(OutputControlSequence);
		mxFree(OutputControlSequence);
	}
}

void IExtInterface::takeInitialStateFromMatlabStruct(
	IExtInterface::SingleStateStruct & IExtInitialStateStruct, 
	mxArray * IExtMatlabInitState, 
	InputArgs & SimulationInputArgs)
{
	int N = SimulationInputArgs.a.size();
	
	// Initializing Irand
	getInputfromStruct<float>(IExtMatlabInitState, "InitialState.Iext.Irand", IExtInitialStateStruct.Irand, 1, "required_size", N);
	
	// Initializing GenState
	{
		bool isNotSingleSeed =
			getInputfromStruct<uint32_t>(IExtMatlabInitState, "InitialState.Iext.GenState", IExtInitialStateStruct.GenState,
				2, "required_size", 1, "no_except");
		if (isNotSingleSeed)
			getInputfromStruct<uint32_t>(IExtMatlabInitState, "InitialState.Iext.GenState", IExtInitialStateStruct.GenState,
				1, "required_size", 4);
	}
}

void IExtInterface::initInternalVariables(
	IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	IExtInterface::InputVarsStruct    & IExtInputVarsStruct, 
	IExtInterface::SingleStateStruct  & IExtInitialStateStruct, 
	InputArgs                         & SimulationInputArgs)
{
	// Aliasing input argument structs
	auto & IntVars = IExtInternalVarsStruct;
	auto & InputVars = IExtInputVarsStruct;
	auto & InitState = IExtInitialStateStruct;
	auto & SimInputArgs = SimulationInputArgs;

	// Initializing Constants
	IntVars.alpha         = InputVars.alpha;
	IntVars.StdDev        = InputVars.StdDev;
	IntVars.OutputControl = InputVars.OutputControl;
	int N = SimInputArgs.a.size();
	int nSteps = SimInputArgs.NoOfms * SimInputArgs.onemsbyTstep;

	// Initializing Random Current Generator from GenState
	XorShiftPlus RandCurrGen;
	XorShiftPlus::StateStruct RandCurrGenState;

	if (InitState.GenState.size() == 1) {
		RandCurrGen = XorShiftPlus(InitState.GenState[0]);
	}
	else if (InitState.GenState.size() == 4) {
		RandCurrGenState.ConvertVecttoState(InitState.GenState);
		RandCurrGen.setstate(RandCurrGenState);
	}

	// Initializing Irand 
	if (InitState.Irand.istrulyempty()) {
		IntVars.Irand.resize(N);
	}
	else {
		IntVars.Irand.assign(InitState.Irand);
	}
	IntVars.Irand.configure(RandCurrGen, XorShiftPlus(), IntVars.alpha);

	// Initializing RandMat and GenMat from initial conditions
	IntVars.RandMat.resize(8192, N);
	IntVars.GenMat.resize(8192, 4);
	
	int LoopLimit = (nSteps >= 8192) ? 8192 : nSteps + 1;
	/* The following Initialization follows as below:

	   In the zeroth Iteration (representative of the 
	   last iteration of previous simulation), RandMat
	   is initialized such that RandMat[0] = Irand at
	   the end of last previous simulation iteration.
	   This is because of the policy that whatever Irand
	   is used at a particular iteration must be avai-
	   lable in RandMat at the end of that iteration.
	   This otherwise translates into the Generate and 
	   then Use policy (as opposed to Use and then Gen-
	   erate)
	   
	   */
	IntVars.RandMat[0] = IntVars.Irand;
	IntVars.Irand.generator1().getstate().ConvertStatetoVect(IntVars.GenMat[0]);
	for (int j = 1; j < LoopLimit; ++j) {
		IntVars.Irand.generate();
		IntVars.RandMat[j] = IntVars.Irand;
		IntVars.Irand.generator1().getstate().ConvertStatetoVect(IntVars.GenMat[j]);
	}

	// Initializing Iext, IextWORand
	IntVars.IextWORand = MexVector<float>(N, 0.0f);
	IntVars.Iext       = MexVector<float>(N, 0.0f);
}

////////////////////////////////////////////////////////
// Output Struct initialization functions 
////////////////////////////////////////////////////////
void IExtInterface::SingleStateStruct::initialize(
	const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	const InternalVars                      & SimulationInternalVars)
{
	GenState = MexVector<uint32_t>(4);
	Irand = MexVector<float>(SimulationInternalVars.N);
}

void IExtInterface::StateOutStruct::initialize(
	const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	const InternalVars                      & SimulationInternalVars)
{
	// Aliasing above funtion parameter structs
	auto &IntVars = IExtInternalVarsStruct;
	auto &SimIntVars = SimulationInternalVars;

	// Aliasing some Simulation Vatiables
	auto & StorageStepSize = SimIntVars.StorageStepSize;
	auto & beta = SimIntVars.beta;
	auto & onemsbyTstep = SimIntVars.onemsbyTstep;
	auto & N = SimIntVars.N;
	auto   nSteps = SimIntVars.NoOfms * SimIntVars.onemsbyTstep;

	// Initializing TimeDimLen
	size_t TimeDimLen;
	if (StorageStepSize) {
		TimeDimLen = (nSteps >= beta) ? (nSteps - beta) / (StorageStepSize*onemsbyTstep) + 1 : 0;	//No. of times (StorageStepSize * onemsbyTstep)|time happens
	}
	else {
		TimeDimLen = nSteps;
	}

	if (IntVars.OutputControl & IExtInterface::OutOps::GEN_STATE_REQ) {
		this->GenStateOut = MexMatrix<uint32_t>(TimeDimLen, 4);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_RAND_REQ) {
		this->IrandOut = MexMatrix<float>(TimeDimLen, N);
	}
}

void IExtInterface::OutputVarsStruct::initialize(
	const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	const InternalVars                      & SimulationInternalVars)
{
	// Aliasing above funtion parameter structs
	auto &IntVars = IExtInternalVarsStruct;
	auto &SimIntVars = SimulationInternalVars;

	// Aliasing some Simulation Vatiables
	auto & StorageStepSize = SimIntVars.StorageStepSize;
	auto & beta            = SimIntVars.beta;
	auto & onemsbyTstep    = SimIntVars.onemsbyTstep;
	auto & N               = SimIntVars.N;
	auto   nSteps          = SimIntVars.NoOfms * SimIntVars.onemsbyTstep;

	// Initializing TimeDimLen
	size_t TimeDimLen;
	if (StorageStepSize) {
		TimeDimLen = (nSteps >= beta) ? (nSteps - beta) / (StorageStepSize*onemsbyTstep) + 1 : 0;	//No. of times (StorageStepSize * onemsbyTstep)|time happens
	}
	else {
		TimeDimLen = nSteps;
	}

	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_WO_RAND_REQ) {
		this->IextWORand = MexMatrix<float>(TimeDimLen, N);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_TOTAL_REQ) {
		this->IextTotal = MexMatrix<float>(TimeDimLen, N);
	}
}

////////////////////////////////////////////////////////
// C++ Output Functions 
////////////////////////////////////////////////////////
void IExtInterface::doSparseOutput(
	IExtInterface::StateOutStruct     & IExtStateOutStruct,
	IExtInterface::OutputVarsStruct   & IExtOutputVarsStruct,
	IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	InternalVars                      & SimulationInternalVars)
{
	// Aliasing above function parameter structs
	auto &IntVars    = IExtInternalVarsStruct;
	auto &OutVars    = IExtOutputVarsStruct;
	auto &StateOut   = IExtStateOutStruct;
	auto &SimIntVars = SimulationInternalVars;

	// Aliasing some simulation variables
	auto &i               = SimIntVars.i;
	auto &beta            = SimIntVars.beta;
	auto &onemsbyTstep    = SimIntVars.onemsbyTstep;
	auto &StorageStepSize = SimIntVars.StorageStepSize;
	auto &DelayRange      = SimIntVars.DelayRange;

	// Initializing some relevant constants
	size_t CurrentInsertPos = (i - beta) / (onemsbyTstep * StorageStepSize);
	size_t IRandIter        = i % 8192;
	size_t QueueSize        = onemsbyTstep * DelayRange;

	// Aliasing some IExtInterface variables
	auto &OutputControl = IntVars.OutputControl;

	// ------------------ OUTPUTTING STATE VARIABLES ------------------ //

	if (IntVars.OutputControl & IExtInterface::OutOps::GEN_STATE_REQ) {
		StateOut.GenStateOut[CurrentInsertPos] = IntVars.GenMat[IRandIter];
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_RAND_REQ) {
		StateOut.IrandOut[CurrentInsertPos] = IntVars.RandMat[IRandIter];
	}

	// ------------------ OUTPUTTING OUTPUT VARIABLES ------------------ //

	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_WO_RAND_REQ) {
		OutVars.IextWORand[CurrentInsertPos] = IntVars.IextWORand;
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_TOTAL_REQ) {
		OutVars.IextTotal[CurrentInsertPos] = IntVars.Iext;
	}
}

void IExtInterface::doFullOutput(
	IExtInterface::StateOutStruct     & IExtStateOutStruct, 
	IExtInterface::OutputVarsStruct   & IExtOutputVarsStruct,
	IExtInterface::InternalVarsStruct & IExtInternalVarsStruct, 
	InternalVars                      & SimulationInternalVars)
{
	// Aliasing above function parameter structs
	auto &IntVars    = IExtInternalVarsStruct;
	auto &OutVars    = IExtOutputVarsStruct;
	auto &StateOut   = IExtStateOutStruct;
	auto &SimIntVars = SimulationInternalVars;

	// Aliasing some simulation variables
	auto &i               = SimIntVars.i;
	auto &beta            = SimIntVars.beta;
	auto &StorageStepSize = SimIntVars.StorageStepSize;

	// Initializing some relevant constants
	size_t CurrentInsertPos = i - 1;
	size_t IRandIter        = i % 8192;

	// Aliasing some IExtInterface variables
	auto &OutputControl = IntVars.OutputControl;

	// ------------------ OUTPUTTING STATE VARIABLES ------------------ //

	if (IntVars.OutputControl & IExtInterface::OutOps::GEN_STATE_REQ) {
		StateOut.GenStateOut[CurrentInsertPos] = IntVars.GenMat[IRandIter];
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_RAND_REQ) {
		StateOut.IrandOut[CurrentInsertPos] = IntVars.RandMat[IRandIter];
	}

	// ------------------ OUTPUTTING OUTPUT VARIABLES ------------------ //

	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_WO_RAND_REQ) {
		OutVars.IextWORand[CurrentInsertPos] = IntVars.IextWORand;
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_TOTAL_REQ) {
		OutVars.IextTotal[CurrentInsertPos] = IntVars.Iext;
	}
}

void IExtInterface::doSingleStateOutput(
	IExtInterface::SingleStateStruct  & IExtSingleStateStruct, 
	IExtInterface::InternalVarsStruct & IExtInternalVarsStruct, 
	InternalVars                      & SimulationInternalVars)
{
	// Aliasing above function parameter structs
	auto &SingleState = IExtSingleStateStruct;
	auto &IntVars     = IExtInternalVarsStruct;
	auto &SimIntVars  = SimulationInternalVars;

	// Aliasing some simulation variables
	auto &i = SimIntVars.i;

	// Initializing some relevant constants
	size_t IRandIter = i % 8192;

	// ------------------ OUTPUTTING STATE VARIABLES ------------------ //

	SingleState.GenState = IntVars.GenMat[IRandIter];
	SingleState.Irand    = IntVars.RandMat[IRandIter];
}

void IExtInterface::doInputVarsOutput(
	IExtInterface::InputVarsStruct    & IExtInInputVarsStruct, 
	IExtInterface::InternalVarsStruct & IExtInternalVarsStruct, 
	InternalVars                      & SimulationInternalVars)
{
	// Aliasing above function parameter structs
	auto & InputVars  = IExtInInputVarsStruct;
	auto & IntVars    = IExtInternalVarsStruct;
	auto & SimIntVars = SimulationInternalVars;

	// Assigning the Input Variables
	InputVars.alpha = IntVars.alpha;
	InputVars.StdDev = IntVars.StdDev;
	// Note that OutputControl for IExtInterface is not so much an input variable
	// as an intermediate variable calculated from an input to the original Simu-
	// lation. Thus, this variablewill not be  returned or passed as input to the
	// Simulation function
}

////////////////////////////////////////////////////////
// Functions performing output to MATLAB Array
////////////////////////////////////////////////////////
mxArrayPtr IExtInterface::putSingleStatetoMATLABStruct(IExtInterface::SingleStateStruct & IExtSingleStateStruct)
{
	const char *FieldNames[] = {
		"GenState",
		"Irand",
		nullptr
	};
	
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & SingleState = IExtSingleStateStruct;

	// Performing output of Single State variables
	mxSetField(ReturnPointer, 0, "GenState", assignmxArray(SingleState.GenState, mxUINT32_CLASS));
	mxSetField(ReturnPointer, 0, "Irand"   , assignmxArray(SingleState.Irand   , mxSINGLE_CLASS));

	return ReturnPointer;
}

mxArrayPtr IExtInterface::putInputVarstoMATLABStruct(IExtInterface::InputVarsStruct & IExtInputVarsStruct)
{
	const char *FieldNames[] = {
		"alpha",
		"StdDev",
		nullptr
	};
	
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & InputVars = IExtInputVarsStruct;
	
	// Performing output of Input variables
	mxSetField(ReturnPointer, 0, "alpha" , assignmxArray(InputVars.alpha , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "StdDev", assignmxArray(InputVars.StdDev, mxSINGLE_CLASS));

	return ReturnPointer;
}

mxArrayPtr IExtInterface::putStateVarstoMATLABStruct(IExtInterface::StateOutStruct & IExtStateOutStruct)
{
	const char *FieldNames[] = {
		"GenState",
		"Irand",
		nullptr
	};

	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & StateVars = IExtStateOutStruct;
	
	// Performing output of Input variables
	mxSetField(ReturnPointer, 0, "GenState", assignmxArray(StateVars.GenStateOut, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "Irand"   , assignmxArray(StateVars.IrandOut   , mxSINGLE_CLASS));

	return ReturnPointer;
}

mxArrayPtr IExtInterface::putOutputVarstoMATLABStruct(IExtInterface::OutputVarsStruct & IExtOutputVarsStruct)
{
	const char *FieldNames[] = {
		"IextTotal",
		"IextWORand",
		nullptr
	};

	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & OutputVars = IExtOutputVarsStruct;

	// Performing output of Input variables
	mxSetField(ReturnPointer, 0, "IextTotal" , assignmxArray(OutputVars.IextTotal , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "IextWORand", assignmxArray(OutputVars.IextWORand, mxSINGLE_CLASS));

	return ReturnPointer;
}

////////////////////////////////////////////////////////
// Iext calculation and update stage 
////////////////////////////////////////////////////////
void IExtInterface::updateIExt(
	IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	InternalVars                      & SimulationInternalVars)
{
	// Aliasing above function parameter structs
	auto &IntVars = IExtInternalVarsStruct;
	auto &SimIntVars = SimulationInternalVars;

	// Initializing Constants
	auto &time         = SimIntVars.Time;
	auto &i            = SimIntVars.i;
	auto &onemsbyTstep = SimIntVars.onemsbyTstep;
	auto nSteps        = SimIntVars.NoOfms * SimIntVars.onemsbyTstep;
	auto &N            = SimIntVars.N;

	// Generating Irand (in the form of RandMat)
	if (i % 8192 == 0) {
		size_t LoopLimit = (i + 8192 <= nSteps) ? 8192 : nSteps - i + 1;
		// The logic for the above is that if i + 8192 > nSteps, then, 
		// this generation will not happen again. In that case, RandMat
		// generation should be done for current iterations plus iter-
		// ations uptil nSteps, hence nSteps - i + 1

		for (int j = 0; j < LoopLimit; ++j) {
			IntVars.Irand.generate();
			IntVars.RandMat[j] = IntVars.Irand;
			IntVars.Irand.generator1().getstate().ConvertStatetoVect(IntVars.GenMat[j]);
		}
	}

	// Generating Iext without Irand using Irand above.
	if (time <= 0.1 * 1000 * onemsbyTstep) {
		for (int j = 0; j < 100 * N / 2000; ++j)
			IntVars.IextWORand[j] = 9.0f;
	}
	else if (time - 0.8 <= 0.015) {	
		for (int j = 0; j < 100 * N / 2000; ++j)
			IntVars.IextWORand[j] = 9.0f;
	}
	else {
		for (int j = 0; j < 100 * N / 2000; ++j)
			IntVars.IextWORand[j] = 9.0f;
	}

	// Generating Total Iext
	for (int j = 0; j < N; ++j) {
		IntVars.Iext[j] = IntVars.IextWORand[j] + IntVars.StdDev*IntVars.RandMat[i % 8192][j];
	}
}
