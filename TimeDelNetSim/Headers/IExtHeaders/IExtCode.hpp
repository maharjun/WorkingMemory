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
			I_RAND_REQ        = (1 << 0),
			GEN_STATE_REQ     = (1 << 1),
			I_EXT_WO_RAND_REQ = (1 << 3),
			I_EXT_TOTAL_REQ   = (1 << 4)
		};
	};

	struct SingleStateStruct {
		MexVector<float>    Irand;
		MexVector<uint32_t> GenState;

		SingleStateStruct() :
			Irand(),
			GenState() {}

		void initialize(
			const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
			const InternalVars                      & SimulationInternalVars);
	};

	struct InputVarsStruct {
		// These are currently the Input Variables
		// State Variables area part of SingleStateStruct

		// Optional Simulation Algorithm Parameters
		float alpha;
		float StdDev;

		size_t OutputControl;

		InputVarsStruct() {}
	};

	struct StateOutStruct {
		MexMatrix<float> IrandOut;
		MexMatrix<uint32_t> GenStateOut;

		StateOutStruct() :
			IrandOut(),
			GenStateOut() {}

		void initialize(
			const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
			const InternalVars                      & SimulationInternalVars);
	};

	struct OutputVarsStruct {
		MexMatrix<float> IextWORand;
		MexMatrix<float> IextTotal;

		OutputVarsStruct() :
			IextWORand(),
			IextTotal() {}

		void initialize(
			const IExtInterface::InternalVarsStruct & IExtInternalVarsStruct,
			const InternalVars                      & SimulationInternalVars);
	};

	struct InternalVarsStruct {

		float alpha;
		float StdDev;

		size_t OutputControl;

		MexMatrix<float> RandMat;
		MexMatrix<uint32_t> GenMat;
		BandLimGaussVect Irand;
		MexVector<float> IextWORand;
		MexVector<float> Iext;

		InternalVarsStruct() :
			RandMat(),
			GenMat(),
			Irand(),
			Iext() {}
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
