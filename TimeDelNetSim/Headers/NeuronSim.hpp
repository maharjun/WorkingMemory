#ifndef NEURONSIM_HPP
#define NEURONSIM_HPP
#include "Network.hpp"
#include "..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\..\MexMemoryInterfacing\Headers\GenericMexIO.hpp"
#include "..\..\MexMemoryInterfacing\Headers\LambdaToFunction.hpp"
#include "..\..\RandomNumGen\Headers\FiltRandomTBB.hpp"
#include ".\IExtHeaders\IExtCode.hpp"


#include <xutility>
#include <stdint.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <tbb\atomic.h>
#include <tbb\parallel_for.h>

#include <emmintrin.h>
#include <smmintrin.h>

#define DEFAULT_STORAGE_STEP 500
#define DEFAULT_STATUS_DISPLAY_STEP 400
using namespace std;

struct OutOps{
	enum {
        CURRENT_QINDS_REQ          = (1 << 1 ),
        FINAL_STATE_REQ            = (1 << 2 ),
        INITIAL_STATE_REQ          = (1 << 3 ),
        I_IN_REQ                   = (1 << 4 ),
        I_TOT_REQ                  = (1 << 5 ),
        LASTSPIKED_NEU_REQ         = (1 << 6 ),
        LASTSPIKED_SYN_REQ         = (1 << 7 ),
        SPIKE_LIST_REQ             = (1 << 8 ),
        SPIKE_QUEUE_REQ            = (1 << 9 ),
        TIME_REQ                   = (1 << 10),
        U_REQ                      = (1 << 11),
        V_REQ                      = (1 << 12),
        WEIGHT_DERIV_REQ           = (1 << 13),
        WEIGHT_REQ                 = (1 << 14),
        NMDA_ACTIVATION_REQ        = (1 << 15),
        NMDA_SYN_EFFECTIVENESS_REQ = (1 << 16),
	};
};

typedef vector<tbb::atomic<long long>, tbb::zero_allocator<long long> > atomicLongVect;
typedef vector<tbb::atomic<int>, tbb::zero_allocator<int> > atomicIntVect;

// Incomplete declarations
struct InputArgs;
struct StateVarsOutStruct;
struct SingleStateStruct;
struct OutputVarsStruct;
struct InternalVars;

struct CurrentUpdate
{
	InternalVars &IntVars;

	CurrentUpdate(InternalVars &IntVars_) :
		IntVars(IntVars_){};
	void operator() (const tbb::blocked_range<int*> &BlockedRange) const;
};
struct NeuronSimulate{
	InternalVars &IntVars;

	NeuronSimulate(
		InternalVars &IntVars_
		) :
		IntVars(IntVars_)
	{};
	void operator() (tbb::blocked_range<int> &Range) const;
};

struct CurrentAttenuate{
	InternalVars &IntVars;

	CurrentAttenuate(
		InternalVars &IntVars_) :
		IntVars(IntVars_){}

	void operator() (tbb::blocked_range<int> &Range) const; 
};

struct SingleStateStruct{
	MexVector<float> Weight;
	MexVector<float> V;
	MexVector<float> U;
	MexVector<float> Iin;
	MexVector<float> WeightDeriv;
	
	// IextInterface state variable component
	IExtInterface::SingleStateStruct IextInterface;

	MexVector<MexVector<int > > SpikeQueue;
	MexVector<int> LSTNeuron;
	MexVector<int> LSTSyn;
	int Time;
	int CurrentQIndex;

	// NMDA Related State Variables
	MexVector<float> NMDA_Activation;
	MexVector<int>   NMDA_SynEffectiveness;

	SingleStateStruct() :
		Weight(),
		V(),
		U(),
		Iin(),
		WeightDeriv(),
		IextInterface(),
		SpikeQueue(),
		LSTNeuron(),
		LSTSyn(),
		NMDA_Activation(),
		NMDA_SynEffectiveness() {}

	void initialize(const InternalVars &);
};

struct InputArgs{
	MexVector<int>   NStart;
	MexVector<int>   NEnd;
	// MexVector<float> Weight; This is a state variable and is initialised 
	//                          and treated as such
	MexVector<float> Delay;
	
	MexVector<float> a;
	MexVector<float> b;
	MexVector<float> c;
	MexVector<float> d;

	MexVector<int> InterestingSyns;

	// IExt Interface Input Variables
	IExtInterface::InputVarsStruct IextInterface;

	// Initial State
	SingleStateStruct InitialState;

	// Compulsory Simulation Parameters
	size_t onemsbyTstep;
	size_t NoOfms;
	size_t DelayRange;

	// Optional Simulation Parameters
	size_t OutputControl;
	size_t StorageStepSize;
	size_t StatusDisplayInterval;
	MexVector<char> OutputControlString;

	// Optional Simulation Algorithm Parameters
	float I0;
	float CurrentDecayFactor;	
	float STDPDecayFactor;
	int   STDPMaxWinLen;
	float MaxSynWeight;
	float W0;

	float NMDA_Ton;             // On  threshold for NMDA Spike
	float NMDA_Toff;            // Off threshold for NMDA Spike
	float NMDA_ActivationDecay; // NMDA Activation decay each ms
	float NMDA_ActivationInc;   // NMDA Activation increment on spike arrival

	InputArgs() :
		NStart(),
		NEnd(),
		Delay(),
		a(),
		b(),
		c(),
		d(),
		InterestingSyns(),
		IextInterface(),
		OutputControlString(),
		InitialState() {}
};

struct InternalVars{
	size_t N;
	size_t NExc;
	size_t M;
	size_t MExc;
	size_t i;       //This is the most important loop index that is definitely a state variable
	                // and plays a crucial role in deciding the index into which the output must be performed
	size_t Time;    // must be initialized befor beta
	size_t beta;    // This is another parameter that plays a rucial role when storing sparsely.
	                // It is the first value of i for which the sparse storage must be done.
	                // goes from 1 to StorageStepSize * onemsbyTstep
	
	// Compulsory Simulation Parameters
	size_t onemsbyTstep;
	size_t NoOfms;
	size_t DelayRange;
	
	// Optional Simulation Parameters
	MexVector<char> OutputControlString;
	size_t OutputControl;
	size_t StorageStepSize;
	const size_t StatusDisplayInterval;

	// Optional Simulation Algorithm Parameters
	const float I0;
	const float CurrentDecayFactor;
	const float STDPDecayFactor;
	const int STDPMaxWinLen;
	const float MaxSynWeight;
	const float W0;
	const float NMDA_Ton;
	const float NMDA_Toff;
	const float NMDA_ActivationDecay;
	const float NMDA_ActivationInc;

	// Scalar State Variables
	// Time is defined earlier for reasons of initialization sequence
	size_t CurrentQIndex;

	// Parameters that control C=Spike Storage Buffering
	size_t CacheBuffering;

	MexVector<Synapse> Network;
	MexVector<Neuron> Neurons;
	MexVector<int> &InterestingSyns;
	MexVector<float> &V;
	MexVector<float> &U;
	atomicLongVect Iin;
	MexVector<float> &WeightDeriv;
	MexVector<float> &NMDA_Activation;
	MexVector<int> &NMDA_SynEffectiveness;

	IExtInterface::InternalVarsStruct IextInterface;

	// Change: remove this. 	
	XorShiftPlus IExtGen;
	size_t CurrentGenNeuron;
	MexVector<MexVector<int> > &SpikeQueue;
	MexVector<int> &LSTNeuron;
	MexVector<int> &LSTSyn;

	MexVector<size_t> AuxArray;						    // Auxillary Array that is an indirection between Network
													    // and an array sorted lexicographically by (NEnd, NStart)
	MexVector<size_t> PreSynNeuronSectionBeg;	        // PreSynNeuronSectionBeg[j] Maintains the list of the 
														// index of the first synapse in Network with NStart = j+1
	MexVector<size_t> PreSynNeuronSectionEnd;	        // PostSynNeuronSectionEnd[j] Maintains the list of the 
														// indices one greater than index of the last synapse in 
														// Network with NStart = j+1

	MexVector<size_t> PostSynNeuronSectionBeg;	        // PostSynNeuronSectionBeg[j] Maintains the list of the 
														// index of the first synapse in AuxArray with NEnd = j+1
	MexVector<size_t> PostSynNeuronSectionEnd;	        // PostSynNeuronSectionEnd[j] Maintains the list of the 
														// indices one greater than index of the last synapse in 
														// AuxArray with NEnd = j+1


	MexVector<float> ExpVect;

	// These vectors are instrumental in the Cache aligned 
	// implementation of Spike Storage
	MexVector<__m128i> BinningBuffer;	//each element is 16 bytes
	MexVector<int> BufferInsertIndex;
	MexVector<int> AddressOffset;

	InternalVars(InputArgs &IArgs) :
		N                     (IArgs.a.size()),
		NExc                  (0),
		M                     (IArgs.NStart.size()),
		MExc                  (0),
		i                     (0),
		Time                  (IArgs.InitialState.Time),
		// beta defined conditionally below
		CurrentQIndex         (IArgs.InitialState.CurrentQIndex),
		OutputControl         (IArgs.OutputControl),
		OutputControlString   (IArgs.OutputControlString),
		StorageStepSize       (IArgs.StorageStepSize),
		StatusDisplayInterval (IArgs.StatusDisplayInterval),
		Network               (N),
		Neurons               (M),
		InterestingSyns       (IArgs.InterestingSyns),
		V                     (IArgs.InitialState.V),
		U                     (IArgs.InitialState.U),		
		Iin                   (N), 
		// Iin is defined separately as an atomic vect.
		WeightDeriv           (IArgs.InitialState.WeightDeriv),
		IExtGen               (),
		CurrentGenNeuron      (0),
		SpikeQueue            (IArgs.InitialState.SpikeQueue),
		LSTNeuron             (IArgs.InitialState.LSTNeuron),
		LSTSyn                (IArgs.InitialState.LSTSyn),
		NMDA_Activation       (IArgs.InitialState.NMDA_Activation),
		NMDA_SynEffectiveness (IArgs.InitialState.NMDA_SynEffectiveness),

		AuxArray                (M),
		PreSynNeuronSectionBeg  (N, -1),
		PreSynNeuronSectionEnd  (N, -1),
		PostSynNeuronSectionBeg (N, -1),
		PostSynNeuronSectionEnd (N, -1),
		ExpVect                 (STDPMaxWinLen),

		BinningBuffer     (CacheBuffering*onemsbyTstep*DelayRange / 4),
		BufferInsertIndex (onemsbyTstep*DelayRange, 0),
		AddressOffset     (onemsbyTstep*DelayRange, 0),

		onemsbyTstep           (IArgs.onemsbyTstep),
		NoOfms                 (IArgs.NoOfms),
		DelayRange             (IArgs.DelayRange),
		CacheBuffering         (128),
		I0                     (IArgs.I0),
		STDPMaxWinLen          (IArgs.STDPMaxWinLen),
		CurrentDecayFactor     (IArgs.CurrentDecayFactor),
		STDPDecayFactor        (IArgs.STDPDecayFactor),
		W0                     (IArgs.W0),
		MaxSynWeight           (IArgs.MaxSynWeight),
		NMDA_Ton               (IArgs.NMDA_Ton),
		NMDA_Toff              (IArgs.NMDA_Toff),
		NMDA_ActivationDecay   (IArgs.NMDA_ActivationDecay),
		NMDA_ActivationInc     (IArgs.NMDA_ActivationInc) {
		
		// Setting up Network and Neurons
		Network.resize(M);
		MexTransform(IArgs.NStart             .begin(), IArgs.NStart             .end(), Network.begin(), FFL([ ](Synapse &Syn, int   &NStart)->void{Syn.NStart        = NStart       ; }));
		MexTransform(IArgs.NEnd               .begin(), IArgs.NEnd               .end(), Network.begin(), FFL([ ](Synapse &Syn, int   &NEnd  )->void{Syn.NEnd          = NEnd         ; }));
		MexTransform(IArgs.InitialState.Weight.begin(), IArgs.InitialState.Weight.end(), Network.begin(), FFL([ ](Synapse &Syn, float &Weight)->void{Syn.Weight        = Weight       ; }));
		MexTransform(IArgs.Delay              .begin(), IArgs.Delay              .end(), Network.begin(), FFL([&](Synapse &Syn, float &Delay )->void{Syn.DelayinTsteps = Delay*IArgs.onemsbyTstep + 0.5f; }));

		Neurons.resize(N);
		MexTransform(IArgs.a.begin(), IArgs.a.end(), Neurons.begin(), FFL([](Neuron &Neu, float &a)->void{Neu.a = a; }));
		MexTransform(IArgs.b.begin(), IArgs.b.end(), Neurons.begin(), FFL([](Neuron &Neu, float &b)->void{Neu.b = b; }));
		MexTransform(IArgs.c.begin(), IArgs.c.end(), Neurons.begin(), FFL([](Neuron &Neu, float &c)->void{Neu.c = c; }));
		MexTransform(IArgs.d.begin(), IArgs.d.end(), Neurons.begin(), FFL([](Neuron &Neu, float &d)->void{Neu.d = d; }));

		// Setting value of beta
		if (StorageStepSize)
			beta = (onemsbyTstep * StorageStepSize) - Time % (onemsbyTstep * StorageStepSize);
		else
			beta = 0;

		// Setting Initial Conditions of V and U
		if (U.istrulyempty()){
			U.resize(N);
			for (int j = 0; j<N; ++j)
				U[j] = Neurons[j].b*(Neurons[j].b - 5.0f - sqrt((5.0f - Neurons[j].b)*(5.0f - Neurons[j].b) - 22.4f)) / 0.08f;
		}
		else if (U.size() != N){
			// GIVE ERROR MESSAGE HERE
			return;
		}

		if (V.istrulyempty()){
			V.resize(N);
			for (int j = 0; j<N; ++j){
				V[j] = (Neurons[j].b - 5.0f - sqrt((5.0f - Neurons[j].b)*(5.0f - Neurons[j].b) - 22.4f)) / 0.08f;
			}
		}
		else if (V.size() != N){
			// GIVE ERROR MESSAGE HEREx
			return;
		}

		// Setting Initial Conditions for INTERNAL CURRENT
		if (IArgs.InitialState.Iin.size() == N){
			for (int j = 0; j < N; ++j){
				Iin[j] = (long long int)(IArgs.InitialState.Iin[j] * (1i64 << 32));
			}
		}
		else if (IArgs.InitialState.Iin.size()){
			// GIVE ERROR MESSAGE HERE
			return;
		}
		//else{
		//	Iin is already initialized to zero by tbb::zero_allocator<long long>
		//}

		// Setting Initial Conditions for Weight Derivative
		if (IArgs.InitialState.WeightDeriv.istrulyempty()){
			WeightDeriv.resize(M, 0.0f);
		}
		else if (IArgs.InitialState.WeightDeriv.size() != M){
			// Return Exception
			return;
		}

		// Initializing InternalVars for IextInternal giving it
		//   1. The IExtInterface::InternalVarsStruct in InternalVars
		//   2. The IExtInterface::InputVarsStruct in InputVars
		//   3. The IExtInterface::SingleStateStruct (Initial State) in InputVars
		//   4. InputVars itself to take parameters such as N, onemsbyTstep, 
		//      Noofms etc.
		
		IExtInterface::initInternalVariables(
			IextInterface,
			IArgs.IextInterface,
			IArgs.InitialState.IextInterface,
			IArgs
		);
		
		// Setting Initial Conditions of SpikeQueue
		if (SpikeQueue.istrulyempty()){
			SpikeQueue = MexVector<MexVector<int> >(onemsbyTstep * DelayRange, MexVector<int>());
		}
		else if (SpikeQueue.size() != onemsbyTstep * DelayRange){
			// GIVE ERROR MESSAGE HERE
			return;
		}

		// Setting Initial Conditions for LSTs
		if (LSTNeuron.istrulyempty()){
			LSTNeuron = MexVector<int>(N, -1);
		}
		else if (LSTNeuron.size() != N){
			//GIVE ERROR MESSAGE HERE
			return;
		}
		if (LSTSyn.istrulyempty()){
			LSTSyn = MexVector<int>(M, -1);
		}
		else if (LSTSyn.size() != M){
			//GIVE ERROR MESSAGE HERE
			return;
		}

		// Setting Initial Conditions for NMDA related variables
		if (NMDA_Activation.istrulyempty()) {
			NMDA_Activation.resize(M, 0.0f);
		}
		if (NMDA_SynEffectiveness.istrulyempty()) {
			NMDA_SynEffectiveness.resize(M, (int)1);
		}
	}
	void DoOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars){
		DoFullOutput(StateOut, OutVars);
		if (StorageStepSize && !(Time % (StorageStepSize*onemsbyTstep))){
			DoSparseOutput(StateOut, OutVars);
		}
	}
	void DoSingleStateOutput(SingleStateStruct &FinalStateOut);
	void DoInputStateOutput(InputArgs &InputStateOut);
private:
	void DoSparseOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars);
	void DoFullOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars);
};

struct OutputVarsStruct{
	MexMatrix<float> WeightOut;
	MexMatrix<float> Itot;

	IExtInterface::OutputVarsStruct IextInterface;

	struct SpikeListStruct{
		MexVector<int> SpikeSynInds;
		MexVector<int> TimeRchdStartInds;
		SpikeListStruct() : SpikeSynInds(), TimeRchdStartInds(){}
	} SpikeList;

	OutputVarsStruct() :
		WeightOut(),
		Itot(),
		SpikeList() {}

	void initialize(const InternalVars &);
};

struct StateVarsOutStruct{
	MexMatrix<float> WeightOut;
	MexMatrix<float> VOut;
	MexMatrix<float> UOut;
	MexMatrix<float> IinOut;
	MexMatrix<float> WeightDerivOut;

	// IExt Interface Output State variables
	IExtInterface::StateOutStruct IextInterface;
	
	MexVector<int> TimeOut;

	MexVector<MexVector<MexVector<int> > > SpikeQueueOut;
	MexVector<int> CurrentQIndexOut;
	MexMatrix<int> LSTNeuronOut;
	MexMatrix<int> LSTSynOut;

	// NMDA State Output Vars
	MexMatrix<float> NMDA_ActivationOut;
	MexMatrix<int> NMDA_SynEffectivenessOut;

	StateVarsOutStruct() :
		WeightOut(),
		VOut(),
		UOut(),
		IinOut(),
		WeightDerivOut(),
		IextInterface(),
		TimeOut(),
		SpikeQueueOut(),
		CurrentQIndexOut(),
		LSTNeuronOut(),
		LSTSynOut(),
		NMDA_ActivationOut(),
		NMDA_SynEffectivenessOut() {}

	void initialize(const InternalVars &);
};

void CountingSort(int N, MexVector<Synapse> &Network, MexVector<int> &indirection);

void SimulateParallel(
	InputArgs &&InputArguments,
	OutputVarsStruct &PureOutputs,
	StateVarsOutStruct &StateVarsOutput,
	SingleStateStruct &FinalStateOutput,
	InputArgs &InitalStateOutput);

#endif