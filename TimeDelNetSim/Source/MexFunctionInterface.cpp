#include <mex.h>
#include <matrix.h>
#undef printf

#include <algorithm>
#include <vector>
#include <cstring>
#include <chrono>
#include <type_traits>

#include "..\Headers\NeuronSim.hpp"

#include "..\Headers\IExtHeaders\IExtCode.hpp"

#include "..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\..\MexMemoryInterfacing\Headers\GenericMexIO.hpp"
#include "..\..\MexMemoryInterfacing\Headers\InterruptHandling.hpp"
#include "..\..\MexMemoryInterfacing\Headers\LambdaToFunction.hpp"

using namespace std;

int getOutputControl(char* OutputControlSequence){
	char * SequenceWord;
	char * NextNonDelim = NULL;
	char * Delims = " -,";
	int OutputControl = 0x00000000;
	SequenceWord = strtok_s(OutputControlSequence, Delims, &NextNonDelim);
	bool AddorRemove; // TRUE for ADD
	while (SequenceWord != NULL) {
		AddorRemove = true;
		if (SequenceWord[0] == '/') {
			AddorRemove = false;
			SequenceWord++;
		}
		if (!_strcmpi(SequenceWord, "Initial"))
			OutputControl |= OutOps::INITIAL_STATE_REQ;
		if (AddorRemove && !_strcmpi(SequenceWord, "VCF"))
			OutputControl |= OutOps::V_REQ | OutOps::U_REQ | OutOps::I_TOT_REQ
		                   | OutOps::FINAL_STATE_REQ;
		if (AddorRemove && !_strcmpi(SequenceWord, "VCWF"))
			OutputControl |= OutOps::V_REQ | OutOps::U_REQ | OutOps::I_TOT_REQ 
			               | OutOps::WEIGHT_REQ
			               | OutOps::FINAL_STATE_REQ;
		if (AddorRemove && !_strcmpi(SequenceWord, "FSF"))
			OutputControl |= OutOps::V_REQ | OutOps::U_REQ 
						   | OutOps::I_IN_1_REQ | OutOps::I_IN_2_REQ
						   | OutOps::WEIGHT_DERIV_REQ 
			               | OutOps::WEIGHT_REQ
						   | OutOps::CURRENT_QINDS_REQ
						   | OutOps::SPIKE_QUEUE_REQ
						   | OutOps::LASTSPIKED_NEU_REQ
						   | OutOps::LASTSPIKED_SYN_REQ
						   | OutOps::FINAL_STATE_REQ;
		if (!_strcmpi(SequenceWord, "V"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::V_REQ : 
					 OutputControl & ~(OutOps::V_REQ);
		if (!_strcmpi(SequenceWord, "U"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::U_REQ : 
					 OutputControl & ~(OutOps::U_REQ);
		if (!_strcmpi(SequenceWord, "Iin1"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_IN_1_REQ : 
					 OutputControl & ~(OutOps::I_IN_1_REQ);
		if (!_strcmpi(SequenceWord, "Iin2"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_IN_2_REQ : 
					 OutputControl & ~(OutOps::I_IN_2_REQ);
		if (!_strcmpi(SequenceWord, "WeightDeriv"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::WEIGHT_DERIV_REQ : 
					 OutputControl & ~(OutOps::WEIGHT_DERIV_REQ);
		if (!_strcmpi(SequenceWord, "Weight"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::WEIGHT_REQ : 
					 OutputControl & ~(OutOps::WEIGHT_REQ);
		if (!_strcmpi(SequenceWord, "CurrentQInds"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::CURRENT_QINDS_REQ : 
					 OutputControl & ~(OutOps::CURRENT_QINDS_REQ);
		if (!_strcmpi(SequenceWord, "SpikeQueue"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::SPIKE_QUEUE_REQ : 
					 OutputControl & ~(OutOps::SPIKE_QUEUE_REQ);
		if (!_strcmpi(SequenceWord, "LSTNeuron"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::LASTSPIKED_NEU_REQ : 
					 OutputControl & ~(OutOps::LASTSPIKED_NEU_REQ);
		if (!_strcmpi(SequenceWord, "LSTSyn"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::LASTSPIKED_SYN_REQ : 
					 OutputControl & ~(OutOps::LASTSPIKED_SYN_REQ);
		if (!_strcmpi(SequenceWord, "Iin"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_IN_REQ : 
					 OutputControl & ~(OutOps::I_IN_REQ);
		if (!_strcmpi(SequenceWord, "Itot"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::I_TOT_REQ : 
					 OutputControl & ~(OutOps::I_TOT_REQ);
		if (!_strcmpi(SequenceWord, "SpikeList"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::SPIKE_LIST_REQ : 
					 OutputControl & ~(OutOps::SPIKE_LIST_REQ);
		if (!_strcmpi(SequenceWord, "Final"))
			OutputControl = AddorRemove ? 
			         OutputControl | OutOps::FINAL_STATE_REQ : 
					 OutputControl & ~(OutOps::FINAL_STATE_REQ);
		SequenceWord = strtok_s(NULL, Delims, &NextNonDelim);
	}
	return OutputControl;
}

void takeInputFromMatlabStruct(mxArray* MatlabInputStruct, InputArgs &InputArgList){

	// Initializing N, M Ensuring that "a" and "NStart" Fields are present
	size_t N = mxGetNumberOfElements(getValidStructField(MatlabInputStruct, "a", MexMemInputOps(true)));
	size_t M = mxGetNumberOfElements(getValidStructField(MatlabInputStruct, "NStart", MexMemInputOps(true)));

	// set Cumpulsory Simulation Parameters
	getInputfromStruct<int>(MatlabInputStruct, "onemsbyTstep", InputArgList.onemsbyTstep, 1, "is_required");
	getInputfromStruct<int>(MatlabInputStruct, "NoOfms"      , InputArgList.NoOfms      , 1, "is_required");
	getInputfromStruct<int>(MatlabInputStruct, "DelayRange"  , InputArgList.DelayRange  , 1, "is_required");

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
	int*        genIntPtr[2];       // Generic int Pointers used around the place to access data
	uint32_t*	genUIntPtr[1];		// Generic unsigned int Pointers used around the place to access data (generator bits)
	short *     genCharPtr;         // Generic short Pointer used around the place to access data (delays specifically)
	mxArray *   genmxArrayPtr;      // Generic mxArray Pointer used around the place to access data

	// Initializing neuron specification structure array Neurons
	getInputfromStruct<float>(MatlabInputStruct, "a", InputArgList.a, 2, "required_size", N, "is_required");
	getInputfromStruct<float>(MatlabInputStruct, "b", InputArgList.b, 2, "required_size", N, "is_required");
	getInputfromStruct<float>(MatlabInputStruct, "c", InputArgList.c, 2, "required_size", N, "is_required");
	getInputfromStruct<float>(MatlabInputStruct, "d", InputArgList.d, 2, "required_size", N, "is_required");

	// Initializing network (Synapse) specification structure array Network
	getInputfromStruct<int>  (MatlabInputStruct, "NStart", InputArgList.NStart , 2, "required_size", M, "is_required");
	getInputfromStruct<int>  (MatlabInputStruct, "NEnd"  , InputArgList.NEnd   , 2, "required_size", M, "is_required");
	getInputfromStruct<float>(MatlabInputStruct, "Weight", InputArgList.Weight , 2, "required_size", M, "is_required");
	getInputfromStruct<float>(MatlabInputStruct, "Delay" , InputArgList.Delay  , 2, "required_size", M, "is_required");

	// Setting Values for Optional Simulation Algorithm Parameters
	getInputfromStruct<float>(MatlabInputStruct, "I0"                 , InputArgList.I0)                 ;
	getInputfromStruct<float>(MatlabInputStruct, "CurrentDecayFactor1", InputArgList.CurrentDecayFactor1);
	getInputfromStruct<float>(MatlabInputStruct, "CurrentDecayFactor2", InputArgList.CurrentDecayFactor2);

	// Initializing Time
	getInputfromStruct<int>(MatlabInputStruct, "InitialState.Time", InputArgList.InitialState.Time);

	// Initializing StorageStepSize
	getInputfromStruct<int>(MatlabInputStruct, "StorageStepSize", InputArgList.StorageStepSize);

	// Initializing StatusDisplayInterval
	getInputfromStruct<int>(MatlabInputStruct, "StatusDisplayInterval", InputArgList.StatusDisplayInterval);

	// Initializing InterestingSyns
	getInputfromStruct<int>(MatlabInputStruct, "InterestingSyns", InputArgList.InterestingSyns);

	// Initializing V, U and Iin1, Iin2
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.V"   , InputArgList.InitialState.V, 1, "required_size", N);
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.U"   , InputArgList.InitialState.U, 1, "required_size", N);
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.Iin1", InputArgList.InitialState.Iin1, 1, "required_size", N);
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.Iin2", InputArgList.InitialState.Iin2, 1, "required_size", N);

	// Initializing WeightDeriv
	getInputfromStruct<float>(MatlabInputStruct, "InitialState.WeightDeriv", InputArgList.InitialState.WeightDeriv, 1, "required_size", M);

	// Initializing CurrentQIndex
	getInputfromStruct<int>(MatlabInputStruct, "InitialState.CurrentQIndex", InputArgList.InitialState.CurrentQIndex);

	// Initializing SpikeQueue
	int SpikeQueueSize = InputArgList.onemsbyTstep * InputArgList.DelayRange;
	getInputfromStruct<int>(MatlabInputStruct, "InitialState.SpikeQueue", InputArgList.InitialState.SpikeQueue, 1, "required_size", SpikeQueueSize);

	// Initializing LastSpikedTimeNeuron
	getInputfromStruct<int>(MatlabInputStruct, "InitialState.LSTNeuron", InputArgList.InitialState.LSTNeuron, 1, "required_size", N);

	// Initializing LastSpikedTimeSyn
	getInputfromStruct<int>(MatlabInputStruct, "InitialState.LSTSyn", InputArgList.InitialState.LSTSyn, 1, "required_size", M);

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
		"SpikeList",
		nullptr
	};

	int NFields = 0;
	for (; FieldNames[NFields] != nullptr; ++NFields);
	mwSize StructArraySize[2] = { 1, 1 };

	mxArray * ReturnPointer = mxCreateStructArray_730(2, StructArraySize, NFields, FieldNames);
	
	// Assigning Weightout
	mxSetField(ReturnPointer, 0, "WeightOut", assignmxArray(Output.WeightOut, mxSINGLE_CLASS));
	// Assigning Iin
	mxSetField(ReturnPointer, 0, "Iin", assignmxArray(Output.Iin, mxSINGLE_CLASS));
	// Assigning Itot
	mxSetField(ReturnPointer, 0, "Itot", assignmxArray(Output.Itot, mxSINGLE_CLASS));
	// Assigning Output variables for IExtInterface
	mxArrayPtr IExtOutVarsStruct = IExtInterface::putOutputVarstoMATLABStruct(Output.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext", IExtOutVarsStruct);

	// Assigning SpikeList
	mxArray * SpikeListStructPtr;
		const char *SpikeListFieldNames[] = {
			"SpikeSynInds",
			"TimeRchdStartInds"
		};
		SpikeListStructPtr = mxCreateStructArray(1, StructArraySize, 2, SpikeListFieldNames);
		Output.SpikeList.SpikeSynInds.trim();
		Output.SpikeList.TimeRchdStartInds.trim();
		mxSetField(SpikeListStructPtr, 0, "SpikeSynInds"     , assignmxArray(Output.SpikeList.SpikeSynInds, mxINT32_CLASS));
		mxSetField(SpikeListStructPtr, 0, "TimeRchdStartInds", assignmxArray(Output.SpikeList.TimeRchdStartInds, mxINT32_CLASS));
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
	mxSetField(ReturnPointer, 0, "V"             , assignmxArray(Output.VOut, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "U"             , assignmxArray(Output.UOut, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "Iin1"          , assignmxArray(Output.Iin1Out, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "Iin2"          , assignmxArray(Output.Iin2Out, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "WeightDeriv"   , assignmxArray(Output.WeightDerivOut, mxSINGLE_CLASS));
	
	// Assigning Iext State variables
	mxArrayPtr mxIExtStateVars = IExtInterface::putStateVarstoMATLABStruct(Output.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext"          , mxIExtStateVars);

	// Assigning current time
	mxSetField(ReturnPointer, 0, "Time"          , assignmxArray(Output.TimeOut, mxINT32_CLASS));

	// Assigning Weight
	mxSetField(ReturnPointer, 0, "Weight"        , assignmxArray(Output.WeightOut, mxSINGLE_CLASS));

	// Assigning Spike Queue Related Shiz
	mxSetField(ReturnPointer, 0, "CurrentQIndex" , assignmxArray(Output.CurrentQIndexOut, mxINT32_CLASS));
	// Assigning SpikeQueue

	mxSetField(ReturnPointer, 0, "SpikeQueue"    , assignmxArray(Output.SpikeQueueOut, mxINT32_CLASS));

	// Assigning Last Spiked Time related information
	mxSetField(ReturnPointer, 0, "LSTNeuron"     , assignmxArray(Output.LSTNeuronOut, mxINT32_CLASS));
	mxSetField(ReturnPointer, 0, "LSTSyn"        , assignmxArray(Output.LSTSynOut, mxINT32_CLASS));

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
	mxSetField(ReturnPointer, 0, "V"                 , assignmxArray(SingleStateList.V, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "U"                 , assignmxArray(SingleStateList.U, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "Iin1"              , assignmxArray(SingleStateList.Iin1, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "Iin2"              , assignmxArray(SingleStateList.Iin2, mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "WeightDeriv"       , assignmxArray(SingleStateList.WeightDeriv, mxSINGLE_CLASS));
	
	// Assigning IExt State variables
	mxArrayPtr mxIExtSingleState = IExtInterface::putSingleStatetoMATLABStruct(SingleStateList.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext"              , mxIExtSingleState);

	// Assigning Time
	if (SingleStateList.Time >= 0)
		mxSetField(ReturnPointer, 0, "Time"          , assignmxArray(SingleStateList.Time, mxINT32_CLASS));
	else
		mxSetField(ReturnPointer, 0, "Time"          , mxCreateNumericMatrix(0, 0, mxINT32_CLASS, mxREAL));

	// Assigning WeightOut
	mxSetField(ReturnPointer, 0, "Weight"            , assignmxArray(SingleStateList.Weight, mxSINGLE_CLASS));

	// Assigning Spike Queue Related Shiz
	if (SingleStateList.CurrentQIndex >= 0)
		mxSetField(ReturnPointer, 0, "CurrentQIndex" , assignmxArray(SingleStateList.CurrentQIndex, mxINT32_CLASS));
	else
		mxSetField(ReturnPointer, 0, "CurrentQIndex" , mxCreateNumericMatrix(0, 0, mxINT32_CLASS, mxREAL));
	mxSetField(ReturnPointer, 0, "SpikeQueue"        , assignmxArray(SingleStateList.SpikeQueue, mxINT32_CLASS));

	// Assigning Last Spiked Time related information
	mxSetField(ReturnPointer, 0, "LSTNeuron"         , assignmxArray(SingleStateList.LSTNeuron, mxINT32_CLASS));
	mxSetField(ReturnPointer, 0, "LSTSyn"            , assignmxArray(SingleStateList.LSTSyn, mxINT32_CLASS));

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
	mxSetField(ReturnPointer, 0, "onemsbyTstep"         , assignmxArray(InputStateStruct.onemsbyTstep         , mxINT32_CLASS));
	mxSetField(ReturnPointer, 0, "NoOfms"               , assignmxArray(InputStateStruct.NoOfms               , mxINT32_CLASS));
	mxSetField(ReturnPointer, 0, "DelayRange"           , assignmxArray(InputStateStruct.DelayRange           , mxINT32_CLASS));
	
	// Assigning Input Vectors (Network, Neurons etc.)
	mxSetField(ReturnPointer, 0, "a"                    , assignmxArray(InputStateStruct.a                    , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "b"                    , assignmxArray(InputStateStruct.b                    , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "c"                    , assignmxArray(InputStateStruct.c                    , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "d"                    , assignmxArray(InputStateStruct.d                    , mxSINGLE_CLASS));

	mxSetField(ReturnPointer, 0, "NStart"               , assignmxArray(InputStateStruct.NStart               , mxINT32_CLASS ));
	mxSetField(ReturnPointer, 0, "NEnd"                 , assignmxArray(InputStateStruct.NEnd                 , mxINT32_CLASS ));
	mxSetField(ReturnPointer, 0, "Weight"               , assignmxArray(InputStateStruct.Weight               , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "Delay"                , assignmxArray(InputStateStruct.Delay                , mxSINGLE_CLASS));
	
	mxSetField(ReturnPointer, 0, "InterestingSyns"      , assignmxArray(InputStateStruct.InterestingSyns      , mxINT32_CLASS ));

	// Assigning Optional Simulation Algorithm Parameters
	mxSetField(ReturnPointer, 0, "I0"                   , assignmxArray(InputStateStruct.I0                   , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "CurrentDecayFactor1"  , assignmxArray(InputStateStruct.CurrentDecayFactor1  , mxSINGLE_CLASS));
	mxSetField(ReturnPointer, 0, "CurrentDecayFactor2"  , assignmxArray(InputStateStruct.CurrentDecayFactor2  , mxSINGLE_CLASS));
	
	// Assigning IExtinterface Input Variables
	mxArrayPtr mxIExtInputVars = IExtInterface::putInputVarstoMATLABStruct(InputStateStruct.IextInterface);
	mxSetField(ReturnPointer, 0, "Iext"                 , mxIExtInputVars);
	
	// Assigning Optional Simulation Parameters
	mxSetField(ReturnPointer, 0, "StorageStepSize"      , assignmxArray(InputStateStruct.StorageStepSize      , mxINT32_CLASS ));
	mxSetField(ReturnPointer, 0, "StatusDisplayInterval", assignmxArray(InputStateStruct.StatusDisplayInterval, mxINT32_CLASS ));
	mxSetField(ReturnPointer, 0, "OutputControl"        , mxCreateString(InputStateStruct.OutputControlString.begin()));
	
	// Assigning initial state
	mxSetField(ReturnPointer, 0, "InitialState"         , putSingleStatetoMatlabStruct(InputStateStruct.InitialState));

	return ReturnPointer;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]){
	// NOTE THAT THERE IS NO DATA VALIDATION AS THIS IS EXPECTED TO HAVE 
	// BEEN DONE IN THE MATLAB SIDE OF THE INTERFACE TO THIS MEX FUNCTION

	// Open Memory Usage Account
	size_t MemAccountKey = MemCounter::OpenMemAccount(size_t(3) << 29);

	InputArgs InputArgList;
	takeInputFromMatlabStruct(prhs[0], InputArgList);

	// Declaring Output Vectors
	OutputVarsStruct PureOutput;
	StateVarsOutStruct StateVarsOutput;
	SingleStateStruct FinalStateOutput;
	InputArgs InputStateOutput;

	// Running Simulation Function.
	chrono::system_clock::time_point TStart = chrono::system_clock::now();
	try{
		EnableInterruptHandling();
		SimulateParallel(
			move(InputArgList),
			PureOutput,
			StateVarsOutput,
			FinalStateOutput,
			InputStateOutput);
		DisableInterruptHandling();
	}
	catch(ExOps::ExCodes A){
		if (A == ExOps::EXCEPTION_MEM_FULL){
		#ifdef MEX_LIB
			char OutputString[256];
			sprintf_s(OutputString, 256, "Mem Limit of %lld MB Exceeded\n", (MemCounter::MemUsageLimit) >> 20);
			mexErrMsgIdAndTxt("CppSimException:MemOverFlow", OutputString);
		#elif defined MEX_EXE
			throw A;
		#endif
		}
	}

	chrono::system_clock::time_point TEnd = chrono::system_clock::now();
	WriteOutput("The Time taken = %d milliseconds\n", chrono::duration_cast<chrono::milliseconds>(TEnd - TStart).count());

	mwSize StructArraySize[2] = { 1, 1 };
	
	plhs[0] = putOutputToMatlabStruct(PureOutput);
	plhs[1] = putStateToMatlabStruct(StateVarsOutput);
	plhs[2] = putSingleStatetoMatlabStruct(FinalStateOutput);

	if (nlhs == 4){
		plhs[3] = putInputStatetoMatlabStruct(InputStateOutput);
	}

	// Close Memory Usage Account
	MemCounter::CloseMemAccount(MemAccountKey);
}