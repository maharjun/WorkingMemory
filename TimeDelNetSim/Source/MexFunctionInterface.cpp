#include <mex.h>
#include <matrix.h>
#undef printf

#include <algorithm>
#include <vector>
#include <cstring>
#include <chrono>
#include <type_traits>

#include "../Headers/NeuronSim.hpp"

#include <ExternalInputCurrent/IExtCode.hpp>

#include <MexMemoryInterfacing/Headers/MexMem.hpp>
#include <MexMemoryInterfacing/Headers/GenericMexIO.hpp>
#include <MexMemoryInterfacing/Headers/InterruptHandling.hpp>
#include <MexMemoryInterfacing/Headers/LambdaToFunction.hpp>
#include <MexMemoryInterfacing/Headers/FlatVectTree/FlatVectTree.hpp>

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

int getOutputControl(char* OutputControlSequence){
	char * SequenceWord;
	char * NextNonDelim = NULL;
	const char * Delims = " -,";
	int OutputControl = 0x00000000;
	SequenceWord = STRTOK_FUNC(OutputControlSequence, Delims, &NextNonDelim);
	bool AddorRemove; // TRUE for ADD
	while (SequenceWord != NULL) {
		AddorRemove = true;
		if (SequenceWord[0] == '/') {
			AddorRemove = false;
			SequenceWord++;
		}
		if (!STRCMPI_FUNC(SequenceWord, "Initial"))
			OutputControl |= OutOps::INITIAL_STATE_REQ;
		if (AddorRemove && !STRCMPI_FUNC(SequenceWord, "VCF"))
			OutputControl |= OutOps::V_REQ | OutOps::U_REQ | OutOps::I_TOT_REQ
		                   | OutOps::FINAL_STATE_REQ;
		if (AddorRemove && !STRCMPI_FUNC(SequenceWord, "VCWF"))
			OutputControl |= OutOps::V_REQ | OutOps::U_REQ | OutOps::I_TOT_REQ 
			               | OutOps::WEIGHT_REQ
			               | OutOps::FINAL_STATE_REQ;
		if (AddorRemove && !STRCMPI_FUNC(SequenceWord, "FSF"))
			OutputControl |= OutOps::V_REQ | OutOps::U_REQ 
						   | OutOps::I_IN_1_REQ | OutOps::I_IN_2_REQ
						   | OutOps::WEIGHT_DERIV_REQ 
			               | OutOps::WEIGHT_REQ
						   | OutOps::CURRENT_QINDS_REQ
						   | OutOps::SPIKE_QUEUE_REQ
						   | OutOps::LASTSPIKED_NEU_REQ
						   | OutOps::LASTSPIKED_SYN_REQ
						   | OutOps::FINAL_STATE_REQ;
		if (!STRCMPI_FUNC(SequenceWord, "V"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::V_REQ : 
					 OutputControl & ~(OutOps::V_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "U"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::U_REQ : 
					 OutputControl & ~(OutOps::U_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "Iin1"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_IN_1_REQ : 
					 OutputControl & ~(OutOps::I_IN_1_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "Iin2"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_IN_2_REQ : 
					 OutputControl & ~(OutOps::I_IN_2_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "WeightDeriv"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::WEIGHT_DERIV_REQ : 
					 OutputControl & ~(OutOps::WEIGHT_DERIV_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "Weight"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::WEIGHT_REQ : 
					 OutputControl & ~(OutOps::WEIGHT_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "CurrentQInds"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::CURRENT_QINDS_REQ : 
					 OutputControl & ~(OutOps::CURRENT_QINDS_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "SpikeQueue"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::SPIKE_QUEUE_REQ : 
					 OutputControl & ~(OutOps::SPIKE_QUEUE_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "LSTNeuron"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::LASTSPIKED_NEU_REQ : 
					 OutputControl & ~(OutOps::LASTSPIKED_NEU_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "LSTSyn"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::LASTSPIKED_SYN_REQ : 
					 OutputControl & ~(OutOps::LASTSPIKED_SYN_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "Iin"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_IN_REQ : 
					 OutputControl & ~(OutOps::I_IN_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "Itot"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_TOT_REQ : 
					 OutputControl & ~(OutOps::I_TOT_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "SpikeList"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::SPIKE_LIST_REQ : 
					 OutputControl & ~(OutOps::SPIKE_LIST_REQ);
		if (!STRCMPI_FUNC(SequenceWord, "Final"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::FINAL_STATE_REQ : 
					 OutputControl & ~(OutOps::FINAL_STATE_REQ);
		SequenceWord = STRTOK_FUNC(NULL, Delims, &NextNonDelim);
	}
	return OutputControl;
}

void takeInputFromMatlabStruct(const mxArray* MatlabInputStruct, InputArgs &InputArgList){
	
	// Initializing N, M Ensuring that "a" and "NStart" Fields are present
	size_t N = mxGetNumberOfElements(getValidStructField(MatlabInputStruct, "a", getInputOps(1, "is_required", "is_nonempty")));
	size_t M = mxGetNumberOfElements(getValidStructField(MatlabInputStruct, "NStart", getInputOps(1, "is_required", "is_nonempty")));

	// set Cumpulsory Simulation Parameters
	getInputfromStruct<int32_t>(MatlabInputStruct, "onemsbyTstep", InputArgList.onemsbyTstep, getInputOps(1, "is_required"));
	getInputfromStruct<int32_t>(MatlabInputStruct, "NoOfms"      , InputArgList.NoOfms      , getInputOps(1, "is_required"));
	getInputfromStruct<int32_t>(MatlabInputStruct, "DelayRange"  , InputArgList.DelayRange  , getInputOps(1, "is_required"));

	// set default values of Optional Simulation Parameters
	InputArgList.StorageStepSize = DEFAULT_STORAGE_STEP;
	InputArgList.OutputControl = 0;
	InputArgList.StatusDisplayInterval = DEFAULT_STATUS_DISPLAY_STEP;
	
	// set default values of Optional Simulation Algorithm Parameters
	InputArgList.I0                  = 1.0f;
	InputArgList.CurrentDecayFactor1 = powf(9.0f / 10, 1.0f / InputArgList.onemsbyTstep);
	InputArgList.CurrentDecayFactor2 = powf(9.0f / (10.0f), 1.0f / (4 * InputArgList.onemsbyTstep));

	// set default values for Scalar State Variables
	InputArgList.InitialState.CurrentQIndex = 0;
	InputArgList.InitialState.Time = 0;

	float*      genFloatPtr[4];     // Generic float Pointers used around the place to access data
	int32_t*    genIntPtr[2];       // Generic signed int32_t Pointers used around the place to access data
	uint32_t*	genUIntPtr[1];		// Generic unsigned int32_t Pointers used around the place to access data (generator bits)
	short *     genCharPtr;         // Generic short Pointer used around the place to access data (delays specifically)
	mxArray *   genmxArrayPtr;      // Generic mxArray Pointer used around the place to access data

	// Initializing neuron specification structure array Neurons
	getInputfromStruct<float>(MatlabInputStruct, "a", InputArgList.a, getInputOps(2, "required_size", N, "is_required"));
	getInputfromStruct<float>(MatlabInputStruct, "b", InputArgList.b, getInputOps(2, "required_size", N, "is_required"));
	getInputfromStruct<float>(MatlabInputStruct, "c", InputArgList.c, getInputOps(2, "required_size", N, "is_required"));
	getInputfromStruct<float>(MatlabInputStruct, "d", InputArgList.d, getInputOps(2, "required_size", N, "is_required"));

	// Initializing network (Synapse) specification structure array Network
	getInputfromStruct<int32_t>(MatlabInputStruct, "NStart", InputArgList.NStart , getInputOps(2, "required_size", M, "is_required"));
	getInputfromStruct<int32_t>(MatlabInputStruct, "NEnd"  , InputArgList.NEnd   , getInputOps(2, "required_size", M, "is_required"));
	getInputfromStruct<float>  (MatlabInputStruct, "Weight", InputArgList.Weight , getInputOps(2, "required_size", M, "is_required"));
	getInputfromStruct<float>  (MatlabInputStruct, "Delay" , InputArgList.Delay  , getInputOps(2, "required_size", M, "is_required"));

	// Setting Values for Optional Simulation Algorithm Parameters
	getInputfromStruct<float>(MatlabInputStruct, "I0"                 , InputArgList.I0)                 ;
	getInputfromStruct<float>(MatlabInputStruct, "CurrentDecayFactor1", InputArgList.CurrentDecayFactor1);
	getInputfromStruct<float>(MatlabInputStruct, "CurrentDecayFactor2", InputArgList.CurrentDecayFactor2);

	// Initializing Time
	getInputfromStruct<int32_t>(MatlabInputStruct, "InitialState.Time", InputArgList.InitialState.Time);

	// Initializing StorageStepSize
	getInputfromStruct<int32_t>(MatlabInputStruct, "StorageStepSize", InputArgList.StorageStepSize);

	// Initializing StatusDisplayInterval
	getInputfromStruct<int32_t>(MatlabInputStruct, "StatusDisplayInterval", InputArgList.StatusDisplayInterval);

	// Initializing InterestingSyns
	getInputfromStruct<int32_t>(MatlabInputStruct, "InterestingSyns", InputArgList.InterestingSyns);

	// Initializing V, U and Iin1, Iin2
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.V"   , InputArgList.InitialState.V, getInputOps(1, "required_size", N));
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.U"   , InputArgList.InitialState.U, getInputOps(1, "required_size", N));
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.Iin1", InputArgList.InitialState.Iin1, getInputOps(1, "required_size", N));
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.Iin2", InputArgList.InitialState.Iin2, getInputOps(1, "required_size", N));

	// Initializing WeightDeriv
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.WeightDeriv", InputArgList.InitialState.WeightDeriv, getInputOps(1, "required_size", M));

	// Initializing CurrentQIndex
	getInputfromStruct<int32_t>(MatlabInputStruct, "InitialState.CurrentQIndex", InputArgList.InitialState.CurrentQIndex);

	// Initializing SpikeQueue
	int SpikeQueueSize = InputArgList.onemsbyTstep * InputArgList.DelayRange;
	getInputfromStruct<int32_t>(MatlabInputStruct, "InitialState.SpikeQueue", InputArgList.InitialState.SpikeQueue, 1, getInputOps(1, "required_size", SpikeQueueSize));

	// Initializing LastSpikedTimeNeuron
	getInputfromStruct<int32_t>(MatlabInputStruct, "InitialState.LSTNeuron", InputArgList.InitialState.LSTNeuron, getInputOps(1, "required_size", N));

	// Initializing LastSpikedTimeSyn
	getInputfromStruct<int32_t>(MatlabInputStruct, "InitialState.LSTSyn", InputArgList.InitialState.LSTSyn, getInputOps(1, "required_size", M));

	// Initializing IExtInterface Input Variables
	IExtInterface::takeInputVarsFromMatlabStruct(InputArgList.IextInterface, MatlabInputStruct, InputArgList);

	// Initializing Iext State variables
	IExtInterface::takeInitialStateFromMatlabStruct(InputArgList.InitialState.IextInterface, MatlabInputStruct, InputArgList);

	// Initializing OutputControl
	// Get OutputControlString and OutputControl Word
	genmxArrayPtr = mxGetField(MatlabInputStruct, 0, "OutputControl");
	if (genmxArrayPtr != NULL && !mxIsEmpty(genmxArrayPtr)){
		char * OutputControlSequence = mxArrayToString(genmxArrayPtr);

		int OutputCtrlSeqLen = std::strlen(OutputControlSequence);
		InputArgList.OutputControlString.resize(OutputCtrlSeqLen + 1);
		InputArgList.OutputControlString.copyArray(0, OutputControlSequence, OutputCtrlSeqLen + 1);

		InputArgList.OutputControl = getOutputControl(OutputControlSequence);
		mxFree(OutputControlSequence);
	}
}

mxArray * putOutputToMatlabStruct(OutputVarsStruct &Output){
	const char *FieldNames[] = { 
		"WeightOut",
		"Iin",
		"Itot",
		"Iext",
		"NoOfSpikes",
		"SpikeList",
		nullptr
	};

	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);
	
	// Assigning Weightout
	mxSetField(ReturnPointer, 0, "WeightOut", assignmxArray(Output.WeightOut));
	// Assigning Iin
	mxSetField(ReturnPointer, 0, "Iin", assignmxArray(Output.Iin));
	// Assigning Itot
	mxSetField(ReturnPointer, 0, "Itot", assignmxArray(Output.Itot));
	// Assigning NoOfSpikes
	mxSetField(ReturnPointer, 0, "NoOfSpikes", assignmxArray(Output.NoOfSpikes));
	// Assigning Output variables for IExtInterface
	mxArrayPtr IExtOutVarsStruct = IExtInterface::putOutputVarstoMATLABStruct(Output.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext", IExtOutVarsStruct);

	// Assigning SpikeList
	mxArray * SpikeListStructPtr;
		const char *SpikeListFieldNames[] = {
			"SpikeSynInds",
			"TimeRchdStartInds",
			"TimeRchd"
		};
		SpikeListStructPtr = mxCreateStructArray(1, StructArraySize, 3, SpikeListFieldNames);
		Output.SpikeList.SpikeSynInds.trim();
		Output.SpikeList.TimeRchdStartInds.trim();
		mxSetField(SpikeListStructPtr, 0, "SpikeSynInds"     , assignmxArray(Output.SpikeList.SpikeSynInds));
		mxSetField(SpikeListStructPtr, 0, "TimeRchdStartInds", assignmxArray(Output.SpikeList.TimeRchdStartInds));
		mxSetField(SpikeListStructPtr, 0, "TimeRchd"         , assignmxArray(Output.SpikeList.TimeRchd         ));
	mxSetField(ReturnPointer, 0, "SpikeList", SpikeListStructPtr);

	return ReturnPointer;
}

mxArray * putStateToMatlabStruct(StateVarsOutStruct &Output){
	const char *FieldNames[] = {
		"V",
		"Iin1",
		"Iin2",
		"WeightDeriv",
		"Iext",
		"Time",
		"U",
		"Weight",
		"CurrentQIndex",
		"SpikeQueue",
		"LSTNeuron",
		"LSTSyn",
		nullptr
	};
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Assigning V, U, Iin, WeightDeriv
	mxSetField(ReturnPointer, 0, "V"             , assignmxArray(Output.VOut));
	mxSetField(ReturnPointer, 0, "U"             , assignmxArray(Output.UOut));
	mxSetField(ReturnPointer, 0, "Iin1"          , assignmxArray(Output.Iin1Out));
	mxSetField(ReturnPointer, 0, "Iin2"          , assignmxArray(Output.Iin2Out));
	mxSetField(ReturnPointer, 0, "WeightDeriv"   , assignmxArray(Output.WeightDerivOut));
	
	// Assigning Iext State variables
	mxArrayPtr mxIExtStateVars = IExtInterface::putStateVarstoMATLABStruct(Output.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext"          , mxIExtStateVars);

	// Assigning current time
	mxSetField(ReturnPointer, 0, "Time"          , assignmxArray(Output.TimeOut));

	// Assigning Weight
	mxSetField(ReturnPointer, 0, "Weight"        , assignmxArray(Output.WeightOut));

	// Assigning Spike Queue Related Shiz
	mxSetField(ReturnPointer, 0, "CurrentQIndex" , assignmxArray(Output.CurrentQIndexOut));
	// Assigning SpikeQueue

	mxSetField(ReturnPointer, 0, "SpikeQueue"    , assignmxArray(Output.SpikeQueueOut));

	// Assigning Last Spiked Time related information
	mxSetField(ReturnPointer, 0, "LSTNeuron"     , assignmxArray(Output.LSTNeuronOut));
	mxSetField(ReturnPointer, 0, "LSTSyn"        , assignmxArray(Output.LSTSynOut));

	return ReturnPointer;
}

mxArray * putSingleStatetoMatlabStruct(SingleStateStruct &SingleStateList){
	const char *FieldNames[] = {
		"V",
		"Iin1",
		"Iin2",
		"WeightDeriv",
		"Iext",
		"Time",
		"U",
		"Weight",
		"CurrentQIndex",
		"SpikeQueue",
		"LSTNeuron",
		"LSTSyn",
		nullptr
	};
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Assigning V, U, Iins, WeightDeriv
	mxSetField(ReturnPointer, 0, "V"                 , assignmxArray(SingleStateList.V));
	mxSetField(ReturnPointer, 0, "U"                 , assignmxArray(SingleStateList.U));
	mxSetField(ReturnPointer, 0, "Iin1"              , assignmxArray(SingleStateList.Iin1));
	mxSetField(ReturnPointer, 0, "Iin2"              , assignmxArray(SingleStateList.Iin2));
	mxSetField(ReturnPointer, 0, "WeightDeriv"       , assignmxArray(SingleStateList.WeightDeriv));
	
	// Assigning IExt State variables
	mxArrayPtr mxIExtSingleState = IExtInterface::putSingleStatetoMATLABStruct(SingleStateList.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext"              , mxIExtSingleState);

	// Assigning Time
	if (SingleStateList.Time >= 0)
		mxSetField(ReturnPointer, 0, "Time"          , assignmxArray<int32_t>(SingleStateList.Time));
	else
		mxSetField(ReturnPointer, 0, "Time"          , mxCreateNumericMatrix(0, 0, mxINT32_CLASS, mxREAL));

	// Assigning WeightOut
	mxSetField(ReturnPointer, 0, "Weight"            , assignmxArray(SingleStateList.Weight));

	// Assigning Spike Queue Related Shiz
	if (SingleStateList.CurrentQIndex >= 0)
		mxSetField(ReturnPointer, 0, "CurrentQIndex" , assignmxArray<int32_t>(SingleStateList.CurrentQIndex));
	else
		mxSetField(ReturnPointer, 0, "CurrentQIndex" , mxCreateNumericMatrix(0, 0, mxINT32_CLASS, mxREAL));
	mxSetField(ReturnPointer, 0, "SpikeQueue"        , assignmxArray(SingleStateList.SpikeQueue));

	// Assigning Last Spiked Time related information
	mxSetField(ReturnPointer, 0, "LSTNeuron"         , assignmxArray(SingleStateList.LSTNeuron));
	mxSetField(ReturnPointer, 0, "LSTSyn"            , assignmxArray(SingleStateList.LSTSyn));

	return ReturnPointer;
}

mxArray * putInputStatetoMatlabStruct(InputArgs &InputStateStruct){
	const char *FieldNames[] = {
		"onemsbyTstep"         ,
		"NoOfms"               ,
		"DelayRange"           ,
		"a"                    ,
		"b"                    ,
		"c"                    ,
		"d"                    ,
		"NStart"               ,
		"NEnd"                 ,
		"Weight"               ,
		"Delay"                ,
		"InterestingSyns"      ,
		"I0"                   ,
		"CurrentDecayFactor1"  ,
		"CurrentDecayFactor2"  ,
		"Iext"                 ,
		"StorageStepSize"      ,
		"OutputControl"        ,
		"StatusDisplayInterval",
		"InitialState"         ,
		nullptr
	};
	
	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);

	// Assigning Compulsory Simulation Parameters
	mxSetField(ReturnPointer, 0, "onemsbyTstep"         , assignmxArray<int32_t>(InputStateStruct.onemsbyTstep         ));
	mxSetField(ReturnPointer, 0, "NoOfms"               , assignmxArray<int32_t>(InputStateStruct.NoOfms               ));
	mxSetField(ReturnPointer, 0, "DelayRange"           , assignmxArray<int32_t>(InputStateStruct.DelayRange           ));
	
	// Assigning Input Vectors (Network, Neurons etc.)
	mxSetField(ReturnPointer, 0, "a"                    , assignmxArray(InputStateStruct.a                    ));
	mxSetField(ReturnPointer, 0, "b"                    , assignmxArray(InputStateStruct.b                    ));
	mxSetField(ReturnPointer, 0, "c"                    , assignmxArray(InputStateStruct.c                    ));
	mxSetField(ReturnPointer, 0, "d"                    , assignmxArray(InputStateStruct.d                    ));

	mxSetField(ReturnPointer, 0, "NStart"               , assignmxArray(InputStateStruct.NStart               ));
	mxSetField(ReturnPointer, 0, "NEnd"                 , assignmxArray(InputStateStruct.NEnd                 ));
	mxSetField(ReturnPointer, 0, "Weight"               , assignmxArray(InputStateStruct.Weight               ));
	mxSetField(ReturnPointer, 0, "Delay"                , assignmxArray(InputStateStruct.Delay                ));
	
	mxSetField(ReturnPointer, 0, "InterestingSyns"      , assignmxArray(InputStateStruct.InterestingSyns      ));

	// Assigning Optional Simulation Algorithm Parameters
	mxSetField(ReturnPointer, 0, "I0"                   , assignmxArray<float>(InputStateStruct.I0                   ));
	mxSetField(ReturnPointer, 0, "CurrentDecayFactor1"  , assignmxArray<float>(InputStateStruct.CurrentDecayFactor1  ));
	mxSetField(ReturnPointer, 0, "CurrentDecayFactor2"  , assignmxArray<float>(InputStateStruct.CurrentDecayFactor2  ));
	
	// Assigning IExtinterface Input Variables
	mxArrayPtr mxIExtInputVars = IExtInterface::putInputVarstoMATLABStruct(InputStateStruct.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext"                 , mxIExtInputVars);
	
	// Assigning Optional Simulation Parameters
	mxSetField(ReturnPointer, 0, "StorageStepSize"      , assignmxArray<int32_t>(InputStateStruct.StorageStepSize));
	mxSetField(ReturnPointer, 0, "StatusDisplayInterval", assignmxArray<int32_t>(InputStateStruct.StatusDisplayInterval));
	mxSetField(ReturnPointer, 0, "OutputControl"        , mxCreateString(InputStateStruct.OutputControlString.begin()));
	
	// Assigning initial state
	mxSetField(ReturnPointer, 0, "InitialState"         , putSingleStatetoMatlabStruct(InputStateStruct.InitialState));

	return ReturnPointer;
}

void mexFunctionCpp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// This function is the function that does the actual work of the MEX
	// Function including performing IO from and to mxArrays, and calling the
	// Simulation function.

	InputArgs InputArgList;
	takeInputFromMatlabStruct(prhs[0], InputArgList);

	// Declaring Output Vectors
	OutputVarsStruct PureOutput;
	StateVarsOutStruct StateVarsOutput;
	SingleStateStruct FinalStateOutput;
	InputArgs InputStateOutput;

	// Running Simulation Function.
	chrono::system_clock::time_point TStart = chrono::system_clock::now();
	
	EnableInterruptHandling();
	SimulateParallel(
		move(InputArgList),
		PureOutput,
		StateVarsOutput,
		FinalStateOutput,
		InputStateOutput);
	DisableInterruptHandling();

	chrono::system_clock::time_point TEnd = chrono::system_clock::now();
	WriteOutput("The Time taken = %d milliseconds\n", chrono::duration_cast<chrono::milliseconds>(TEnd - TStart).count());

	mwSize StructArraySize[2] = { 1, 1 };
	
	plhs[0] = putOutputToMatlabStruct(PureOutput);
	plhs[1] = putStateToMatlabStruct(StateVarsOutput);
	plhs[2] = putSingleStatetoMatlabStruct(FinalStateOutput);

	if (nlhs == 4){
		plhs[3] = putInputStatetoMatlabStruct(InputStateOutput);
	}
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