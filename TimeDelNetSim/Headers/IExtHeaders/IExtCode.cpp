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
				IExtInterface::OutOps::I_EXT_GEN_STATE_REQ |
				IExtInterface::OutOps::I_EXT_REQ |
				IExtInterface::OutOps::I_EXT_NEURON_REQ |
				IExtInterface::OutOps::I_RAND_NEURON_REQ;
		}

		if (iEqual(OutputControlOptionParts[0], "IExt") || iEqual(OutputControlOptionParts[0], "/IExt")
			&& OutputControlOptionParts.size() == 2) {
			// Ascertain Add or Remove
			if (OutputControlOptionParts[0][0] == '/') {
				AddorRemove = false;
			}

			// Check out different Output Option cases
			if (iEqual(OutputControlOptionParts[1], "IExtGenState")) { 
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_EXT_GEN_STATE_REQ: 
					OutputControlWord & ~IExtInterface::OutOps::I_EXT_GEN_STATE_REQ;
			}
			if (iEqual(OutputControlOptionParts[1], "Iext")) {
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_EXT_REQ: 
					OutputControlWord & ~IExtInterface::OutOps::I_EXT_REQ;
			}
			if (iEqual(OutputControlOptionParts[1], "IExtNeuron")) {
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_EXT_NEURON_REQ:
					OutputControlWord & ~IExtInterface::OutOps::I_EXT_NEURON_REQ;
			}
			if (iEqual(OutputControlOptionParts[1], "IRandNeuron")) {
				OutputControlWord = (AddorRemove) ?
					OutputControlWord |  IExtInterface::OutOps::I_RAND_NEURON_REQ: 
					OutputControlWord & ~IExtInterface::OutOps::I_RAND_NEURON_REQ;
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
	IExtInputVarsStruct.IRandDecayFactor = 2.0f / 3;
	IExtInputVarsStruct.IRandAmplitude = 20.0f;
	IExtInputVarsStruct.IExtAmplitude  = 20.0f;

	// Taking input for Optional Simulation Algorithm Parameters
	getInputfromStruct<float>(IExtMatlabInputStruct, "Iext.IRandDecayFactor", IExtInputVarsStruct.IRandDecayFactor);
	getInputfromStruct<float>(IExtMatlabInputStruct, "Iext.IRandAmplitude"  , IExtInputVarsStruct.IRandAmplitude  );
	getInputfromStruct<float>(IExtMatlabInputStruct, "Iext.IExtAmplitude"   , IExtInputVarsStruct.IExtAmplitude   );

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
	getInputfromStruct<float>(IExtMatlabInitState, "InitialState.Iext.Iext", IExtInitialStateStruct.Iext, 1, "required_size", N);
	
	// Initializing IExtGenState
	{
		bool isNotSingleSeed =
			getInputfromStruct<uint32_t>(IExtMatlabInitState, "InitialState.Iext.IExtGenState", IExtInitialStateStruct.IExtGenState,
				2, "required_size", 1, "no_except");
		if (isNotSingleSeed)
			getInputfromStruct<uint32_t>(IExtMatlabInitState, "InitialState.Iext.IExtGenState", IExtInitialStateStruct.IExtGenState,
				1, "required_size", 4);
	}

	// Initializing IExtNeuron
	IExtInitialStateStruct.IExtNeuron = 0;
	getInputfromStruct<int>(IExtMatlabInitState, "InitialState.Iext.IExtNeuron", IExtInitialStateStruct.IExtNeuron);

	// Initializing IRandNeuron
	IExtInitialStateStruct.IRandNeuron = 0;
	getInputfromStruct<int>(IExtMatlabInitState, "InitialState.Iext.IRandNeuron", IExtInitialStateStruct.IRandNeuron);

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

	// Initializing Input Vars
	IntVars.IRandDecayFactor = InputVars.IRandDecayFactor;
	IntVars.IRandAmplitude   = InputVars.IRandAmplitude;
	IntVars.IExtAmplitude    = InputVars.IExtAmplitude;
	IntVars.OutputControl    = InputVars.OutputControl;

	// ---------- INITIALIZING STATE VARIABLES ---------- //

	// Initializing required constants
	int N = SimInputArgs.a.size();
	int nSteps = SimInputArgs.NoOfms * SimInputArgs.onemsbyTstep;

	// Initializing IExtGen from IExtGenState
	XorShiftPlus::StateStruct RandCurrGenState;

	if (InitState.IExtGenState.size() == 1) {
		IntVars.IExtGen = XorShiftPlus(InitState.IExtGenState[0]);
	}
	else if (InitState.IExtGenState.size() == 4) {
		RandCurrGenState.ConvertVecttoState(InitState.IExtGenState);
		IntVars.IExtGen.setstate(RandCurrGenState);
	}

	// Initializing Iext
	if (InitState.Iext.istrulyempty())
		IntVars.Iext.resize(N, 0.0f);
	else
		IntVars.Iext = InitState.Iext;
	
	// Initializing IExtNeuron
	IntVars.IExtNeuron = InitState.IExtNeuron;

	// Initializing IRandNeuron
	IntVars.IRandNeuron = InitState.IRandNeuron;

}

////////////////////////////////////////////////////////
// Output Struct initialization functions 
////////////////////////////////////////////////////////
void IExtInterface::SingleStateStruct::initialize(
	const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	const InternalVars                      & SimulationInternalVars)
{
	IExtGenState = MexVector<uint32_t>(4);
	Iext         = MexVector<float>(SimulationInternalVars.N);
	IExtNeuron   = 0;
	IRandNeuron  = 0;
}

void IExtInterface::StateOutStruct::initialize(
	const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
	const InternalVars                      & SimulationInternalVars)
{
	// Aliasing above funtion parameter structs
	auto & IntVars = IExtInternalVarsStruct;
	auto & SimIntVars = SimulationInternalVars;

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

	// Initializing state variables output struct based on output options
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_GEN_STATE_REQ) {
		this->IExtGenStateOut = MexMatrix<uint32_t>(TimeDimLen, 4);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_REQ) {
		this->IextOut = MexMatrix<float>(TimeDimLen, N);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_NEURON_REQ) {
		this->IExtNeuronOut = MexVector<int>(TimeDimLen);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_RAND_NEURON_REQ) {
		this->IRandNeuronOut = MexVector<int>(TimeDimLen);
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

	// Initializing Output Variables according to output options
	// Currently there are no Output variables
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
	size_t QueueSize        = onemsbyTstep * DelayRange;

	// Aliasing some IExtInterface variables
	auto &OutputControl = IntVars.OutputControl;

	// ------------------ OUTPUTTING STATE VARIABLES ------------------ //

	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_GEN_STATE_REQ) {
		IntVars.IExtGen.getstate().ConvertStatetoVect(StateOut.IExtGenStateOut[CurrentInsertPos]);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_REQ) {
		StateOut.IextOut[CurrentInsertPos] = IntVars.Iext;
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_NEURON_REQ) {
		StateOut.IExtNeuronOut[CurrentInsertPos] = IntVars.IExtNeuron;
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_RAND_NEURON_REQ) {
		StateOut.IRandNeuronOut[CurrentInsertPos] = IntVars.IRandNeuron;
	}
	// ------------------ OUTPUTTING OUTPUT VARIABLES ------------------ //

	// No output variables
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

	// Aliasing some IExtInterface variables
	auto &OutputControl = IntVars.OutputControl;

	// ------------------ OUTPUTTING STATE VARIABLES ------------------ //

	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_GEN_STATE_REQ) {
		IntVars.IExtGen.getstate().ConvertStatetoVect(StateOut.IExtGenStateOut[CurrentInsertPos]);
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_REQ) {
		StateOut.IextOut[CurrentInsertPos] = IntVars.Iext;
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_EXT_NEURON_REQ) {
		StateOut.IExtNeuronOut[CurrentInsertPos] = IntVars.IExtNeuron;
	}
	if (IntVars.OutputControl & IExtInterface::OutOps::I_RAND_NEURON_REQ) {
		StateOut.IRandNeuronOut[CurrentInsertPos] = IntVars.IRandNeuron;
	}
	// ------------------ OUTPUTTING OUTPUT VARIABLES ------------------ //

	// No output variables
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

	// ------------------ OUTPUTTING STATE VARIABLES ------------------ //

	IntVars.IExtGen.getstate().ConvertStatetoVect(SingleState.IExtGenState);
	SingleState.Iext        = IntVars.Iext;
	SingleState.IExtNeuron  = IntVars.IExtNeuron;
	SingleState.IRandNeuron = IntVars.IRandNeuron;
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
	InputVars.IRandDecayFactor = IntVars.IRandDecayFactor;
	InputVars.IRandAmplitude   = IntVars.IRandAmplitude;
	InputVars.IExtAmplitude    = IntVars.IExtAmplitude;

	// Note that OutputControl for IExtInterface is not so much an input variable
	// as an intermediate variable calculated from an input to the original Simu-
	// lation. Thus, this variable will not be returned or passed as input to the
	// Simulation function
}

////////////////////////////////////////////////////////
// Functions performing output to MATLAB Array
////////////////////////////////////////////////////////
mxArrayPtr IExtInterface::putSingleStatetoMATLABStruct(IExtInterface::SingleStateStruct & IExtSingleStateStruct)
{
	const char *FieldNames[] = {
		"IExtGenState",
		"Iext"        ,
		"IExtNeuron"  ,
		"IRandNeuron" ,
		nullptr
	};
	
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & SingleState = IExtSingleStateStruct;

	// Performing output of Single State variables
	mxSetField(ReturnPointer, 0, "IExtGenState", assignmxArray(SingleState.IExtGenState, mxUINT32_CLASS));
	mxSetField(ReturnPointer, 0, "Iext"        , assignmxArray(SingleState.Iext        , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "IExtNeuron"  , assignmxArray(SingleState.IExtNeuron  , mxINT32_CLASS));
	mxSetField(ReturnPointer, 0, "IRandNeuron" , assignmxArray(SingleState.IRandNeuron , mxINT32_CLASS));

	return ReturnPointer;
}

mxArrayPtr IExtInterface::putInputVarstoMATLABStruct(IExtInterface::InputVarsStruct & IExtInputVarsStruct)
{
	const char *FieldNames[] = {
		"IRandDecayFactor",
		"IRandAmplitude"  ,
		"IExtAmplitude"   ,
		nullptr
	};
	
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & InputVars = IExtInputVarsStruct;
	
	// Performing output of Input variables
	mxSetField(ReturnPointer, 0, "IRandDecayFactor", assignmxArray(InputVars.IRandDecayFactor, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "IRandAmplitude"  , assignmxArray(InputVars.IRandAmplitude  , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "IExtAmplitude"   , assignmxArray(InputVars.IExtAmplitude   , mxSINGLE_CLASS));

	return ReturnPointer;
}

mxArrayPtr IExtInterface::putStateVarstoMATLABStruct(IExtInterface::StateOutStruct & IExtStateOutStruct)
{
	const char *FieldNames[] = {
		"IExtGenState",
		"Iext"        ,
		"IExtNeuron"  ,
		"IRandNeuron" ,
		nullptr
	};

	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & StateVars = IExtStateOutStruct;
	
	// Performing output of Input variables
	mxSetField(ReturnPointer, 0, "IExtGenState", assignmxArray(StateVars.IExtGenStateOut, mxUINT32_CLASS));
	mxSetField(ReturnPointer, 0, "Iext"        , assignmxArray(StateVars.IextOut        , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "IExtNeuron"  , assignmxArray(StateVars.IExtNeuronOut  , mxINT32_CLASS));
	mxSetField(ReturnPointer, 0, "IRandNeuron" , assignmxArray(StateVars.IRandNeuronOut , mxINT32_CLASS));

	return ReturnPointer;
}

mxArrayPtr IExtInterface::putOutputVarstoMATLABStruct(IExtInterface::OutputVarsStruct & IExtOutputVarsStruct)
{
	const char *FieldNames[] = {
		nullptr
	};

	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Aliasing input parameter structs
	auto & OutputVars = IExtOutputVarsStruct;

	// Performing output of Input variables
	// No output variables

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
	auto &N            = SimIntVars.N;
	auto &Time         = SimIntVars.Time;
	auto &onemsbyTstep = SimIntVars.onemsbyTstep;

	// Resetting IExt
	if (IntVars.IRandNeuron > 0)
		IntVars.Iext[IntVars.IRandNeuron - 1] = 0;
	
	if (IntVars.IExtNeuron > 0) {
		IntVars.Iext[IntVars.IExtNeuron - 1] = 0;
		IntVars.IExtNeuron = 0;
	}
	

	// Random Neuron Selection once every time step (in this case ms)
	// Added to deliberate external current
	//if (Time  < 100*1000*onemsbyTstep)
	if (Time % 10000 < 1000) {
		if (Time % 100 < 60) {
			IntVars.IExtNeuron = Time % 100 + 1;
			IntVars.Iext[IntVars.IExtNeuron - 1] += IntVars.IExtAmplitude;
		}
	}
	IntVars.IRandNeuron = (IntVars.IExtGen() % (int)(3.333f*N)) + 1;
	
	if (IntVars.IRandNeuron <= N) {
		IntVars.Iext[IntVars.IRandNeuron - 1] += IntVars.IRandAmplitude;
	}
	else {
		IntVars.IRandNeuron = 0;
	}

}
