#ifndef I_EXT_CODE_HPP
#define I_EXT_CODE_HPP

#include <matrix.h>
#include <mex.h>
#undef printf

#include "..\..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\..\..\RandomNumGen\Headers\FiltRandomTBB.hpp"

struct InternalVars;
struct InputArgs;

namespace IExtInterface
{
	// Forward declaration of all struct types
	struct SingleStateStruct;
	struct InputVarsStruct;
	struct StateOutStruct;
	struct OutputVarsStruct;
	struct InternalVarsStruct;

	struct OutOps {
		enum {
			I_EXT_REQ           = (1 << 0),
			I_EXT_GEN_STATE_REQ = (1 << 1),
			I_EXT_NEURON_REQ    = (1 << 2),
		};
	};

	struct SingleStateStruct {
		MexVector<float>    Iext;
		MexVector<uint32_t> IExtGenState;
		int                 IExtNeuron;

		SingleStateStruct() :
			Iext(),
			IExtGenState(),
			IExtNeuron() {}

		void initialize(
			const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
			const InternalVars                      & SimulationInternalVars);
	};

	struct InputVarsStruct {
		// These are currently the Input Variables
		// State Variables area part of SingleStateStruct

		// Optional Simulation Algorithm Parameters
		float IExtScaleFactor;
		float IExtDecayFactor;

		size_t OutputControl;

		InputVarsStruct() {}
	};

	struct StateOutStruct {
		MexMatrix<float>    IextOut;
		MexMatrix<uint32_t> IExtGenStateOut;
		MexVector<int>      IExtNeuronOut;

		StateOutStruct() :
			IextOut(),
			IExtGenStateOut(),
			IExtNeuronOut() {}

		void initialize(
			const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
			const InternalVars                      & SimulationInternalVars);
	};

	struct OutputVarsStruct {
		
		OutputVarsStruct() {}

		void initialize(
			const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
			const InternalVars                      & SimulationInternalVars);
	};

	struct InternalVarsStruct {

		float IExtScaleFactor;
		float IExtDecayFactor;

		size_t OutputControl;

		MexVector<float> Iext;
		XorShiftPlus IExtGen;
		size_t IExtNeuron;

		InternalVarsStruct() :
			Iext(),
			IExtGen(), 
			IExtNeuron(0) {}
	};

	size_t getOutputControl(char *OutputControlString);

	// Input and Initialization Functions
	void takeInputVarsFromMatlabStruct(
		IExtInterface::InputVarsStruct &IExtInputVarsStruct,
		mxArray * MatlabInputStruct,
		InputArgs &SimulationInputArgs);

	void takeInitialStateFromMatlabStruct(
		IExtInterface::SingleStateStruct &IExtInitialStateStruct,
		mxArray * MatlabInputStruct,
		InputArgs &SimulationInputArgs);

	void initInternalVariables(
		IExtInterface::InternalVarsStruct &IExtInternalVarsStruct,
		IExtInterface::InputVarsStruct    &IExtInputVarsStruct,
		IExtInterface::SingleStateStruct  &IExtInitialStateStruct,
		InputArgs &SimulationInputArgs);

	// C++ Output Functions
	void doSparseOutput(
		IExtInterface::StateOutStruct     &IExtStateOutStruct,
		IExtInterface::OutputVarsStruct   &IExtOutputVarsStruct,
		IExtInterface::InternalVarsStruct &IExtInternalVarsStruct,
		InternalVars                      &SimulationInternalVars);

	void doFullOutput(
		IExtInterface::StateOutStruct     &IExtStateOutStruct,
		IExtInterface::OutputVarsStruct   &IExtOutputVarsStruct,
		IExtInterface::InternalVarsStruct &IExtInternalVarsStruct,
		InternalVars                      &SimulationInternalVars);

	void doSingleStateOutput(
		IExtInterface::SingleStateStruct  &IExtSingleStateStruct,
		IExtInterface::InternalVarsStruct &IExtInternalVarsStruct,
		InternalVars                      &SimulationInternalVars);

	void doInputVarsOutput(
		IExtInterface::InputVarsStruct    &IExtInInputVarsStruct,
		IExtInterface::InternalVarsStruct &IExtInternalVarsStruct,
		InternalVars                      &SimulationInternalVars);

	// Functions performing output to MATLAB Array
	mxArrayPtr putSingleStatetoMATLABStruct (IExtInterface::SingleStateStruct & IExtSingleStateStruct);
	mxArrayPtr putInputVarstoMATLABStruct   (IExtInterface::InputVarsStruct   & IExtInputVarsStruct);
	mxArrayPtr putStateVarstoMATLABStruct   (IExtInterface::StateOutStruct    & IExtStateOutStruct);
	mxArrayPtr putOutputVarstoMATLABStruct  (IExtInterface::OutputVarsStruct  & IExtOutputVarsStruct);

	// IExt Calculating and updating function
	void updateIExt(
		IExtInterface::InternalVarsStruct &IExtInternalVarsStruct,
		InternalVars                      &SimulationInternalVars);
};

#endif
