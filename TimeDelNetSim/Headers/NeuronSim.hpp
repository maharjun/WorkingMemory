#ifndef NEURONSIM_HPP
#define NEURONSIM_HPP
#include "Network.hpp"
#include "..\..\MexMemoryInterfacing\Headers\MexMem.hpp"
#include "..\..\RandomNumGen\Headers\FiltRandomTBB.hpp"
#include <mex.h>
#include <matrix.h>
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
		V_REQ               = (1 << 0), 
		I_IN_1_REQ          = (1 << 1), 
		I_IN_2_REQ          = (1 << 2), 
		I_IN_REQ            = (1 << 3), 
		I_RAND_REQ          = (1 << 4), 
		WEIGHT_DERIV_REQ    = (1 << 5), 
		GEN_STATE_REQ       = (1 << 6), 
		TIME_REQ            = (1 << 7), 
		U_REQ               = (1 << 8), 
		WEIGHT_REQ          = (1 << 9), 
		CURRENT_QINDS_REQ   = (1 << 10), 
		SPIKE_QUEUE_REQ     = (1 << 11), 
		LASTSPIKED_NEU_REQ  = (1 << 12), 
		LASTSPIKED_SYN_REQ  = (1 << 13), 
		I_TOT_REQ           = (1 << 14), 
		SPIKE_LIST_REQ      = (1 << 15), 
		INITIAL_STATE_REQ   = (1 << 16), 
		FINAL_STATE_REQ     = (1 << 17)
	};
};

typedef vector<tbb::atomic<long long>, tbb::zero_allocator<long long> > atomicLongVect;
typedef vector<tbb::atomic<int>, tbb::zero_allocator<int> > atomicIntVect;

// Incomplete declarations
struct InputArgs;
struct StateVarsOutStruct;
struct SingleStateStruct;
struct FinalStateStruct;
struct InitialStateStruct;
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


struct InputArgs{
	static void IExtFunc(float, MexVector<float> &);
	MexVector<Synapse> Network;
	MexVector<Neuron> Neurons;
	MexVector<int> InterestingSyns;
	MexVector<float> V;
	MexVector<float> U;
	MexVector<float> Iin1;
	MexVector<float> Iin2;
	MexVector<float> WeightDeriv;
	
	MexVector<uint32_t> GenState;
	MexVector<float> Irand;
	
	MexVector<MexVector<int> > SpikeQueue;
	MexVector<int> LSTNeuron;
	MexVector<int> LSTSyn;

	// Compulsory Simulation Parameters
	size_t onemsbyTstep;
	size_t NoOfms;
	size_t DelayRange;

	// Optional Simulation Parameters
	size_t OutputControl;
	size_t StorageStepSize;
	size_t StatusDisplayInterval;

	// Optional Simulation Algorithm Parameters
	float I0;
	float CurrentDecayFactor1, CurrentDecayFactor2;
	float alpha;
	float StdDev;

	// Scalar State Variables
	size_t Time;
	size_t CurrentQIndex;

	InputArgs() :
		Network(),
		Neurons(),
		InterestingSyns(),
		V(),
		U(),
		Iin1(),
		Iin2(),
		WeightDeriv(),
		GenState(),
		Irand(),
		SpikeQueue(),
		LSTNeuron(),
		LSTSyn() {}
};

struct InternalVars{
	size_t N;
	size_t M;
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
	size_t OutputControl;
	size_t StorageStepSize;
	const size_t StatusDisplayInterval;

	// Optional Simulation Algorithm Parameters
	const float I0;
	const float CurrentDecayFactor1, CurrentDecayFactor2;
	const float alpha;
	const float StdDev;

	// Scalar State Variables
	// Time is defined earlier for reasons of initialization sequence
	size_t CurrentQIndex;

	// Parameters that control C=Spike Storage Buffering
	size_t CacheBuffering;

	MexVector<Synapse> &Network;
	MexVector<Neuron> &Neurons;
	MexVector<int> &InterestingSyns;
	MexVector<float> &V;
	MexVector<float> &U;
	atomicLongVect Iin1;
	atomicLongVect Iin2;
	MexVector<float> WeightDeriv;
	BandLimGaussVect Irand;
	MexMatrix<float> RandMat;
	MexMatrix<uint32_t> GenMat;
	MexVector<float> Iext;

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

	// These vectors are instrumental in the Cache aligned 
	// implementation of Spike Storage
	MexVector<__m128i> BinningBuffer;	//each element is 16 bytes
	MexVector<int> BufferInsertIndex;
	MexVector<int> AddressOffset;

	InternalVars(InputArgs &IArgs) :
		N                     (IArgs.Neurons.size()),
		M                     (IArgs.Network.size()),
		i                     (0),
		Time                  (IArgs.Time),
		// beta defined conditionally below
		CurrentQIndex         (IArgs.CurrentQIndex),
		OutputControl         (IArgs.OutputControl),
		StorageStepSize       (IArgs.StorageStepSize),
		StatusDisplayInterval (IArgs.StatusDisplayInterval),
		Network               (IArgs.Network),
		Neurons               (IArgs.Neurons),
		InterestingSyns       (IArgs.InterestingSyns),
		V                     (IArgs.V),
		U                     (IArgs.U),
		Iin1(N),	// Iin is defined separately as an atomic vect.
		Iin2(N),
		WeightDeriv(M),
		Irand(),	// Irand defined separately.
		RandMat(8192, N),
		GenMat(8192, 4),
		Iext(N, 0.0f),
		SpikeQueue            (IArgs.SpikeQueue),
		LSTNeuron             (IArgs.LSTNeuron),
		LSTSyn                (IArgs.LSTSyn),
		AuxArray                (M),
		PreSynNeuronSectionBeg  (N, -1),
		PreSynNeuronSectionEnd  (N, -1),
		PostSynNeuronSectionBeg (N, -1),
		PostSynNeuronSectionEnd (N, -1),

		BinningBuffer     (CacheBuffering*onemsbyTstep*DelayRange / 4),
		BufferInsertIndex (onemsbyTstep*DelayRange, 0),
		AddressOffset     (onemsbyTstep*DelayRange, 0),

		onemsbyTstep       (IArgs.onemsbyTstep),
		NoOfms             (IArgs.NoOfms),
		DelayRange         (IArgs.DelayRange),
		CacheBuffering     (128),
		I0                 (IArgs.I0),
		CurrentDecayFactor1(IArgs.CurrentDecayFactor1),
		CurrentDecayFactor2(IArgs.CurrentDecayFactor2),
		alpha              (IArgs.alpha), 
		StdDev             (IArgs.StdDev)
	{

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

		// Setting Initial Conditions for INTERNAL CURRENT 1
		if (IArgs.Iin1.size() == N){
			for (int j = 0; j < N; ++j){
				Iin1[j] = (long long int)(IArgs.Iin1[j] * (1i64 << 32));
			}
		}
		else if (IArgs.Iin1.size()){
			// GIVE ERROR MESSAGE HERE
			return;
		}
		//else{
		//	Iin1 is already initialized to zero by tbb::zero_allocator<long long>
		//}

		// Setting Initial Conditions for INTERNAL CURRENT 2
		if (IArgs.Iin2.size() == N){
			for (int j = 0; j < N; ++j){
				Iin2[j] = (long long int)(IArgs.Iin2[j] * (1i64 << 32));
		}
		}
		else if (IArgs.Iin2.size()){
			// GIVE ERROR MESSAGE HERE
			return;
		}
		//else{
		//	Iin2 is already initialized to zero by tbb::zero_allocator<long long>
		//}

		// Setting Initial Conditions for Weight Derivative
		if (IArgs.WeightDeriv.istrulyempty()){
			for (auto WeightDerivElem : WeightDeriv) WeightDerivElem = 0;
		}
		else if (IArgs.WeightDeriv.size() == M){
			WeightDeriv = IArgs.WeightDeriv;
		}
		else{
			// Return Exception
			return;
		}

		// Setting up IRand and corresponding Random Generators.
		XorShiftPlus Gen1;
		XorShiftPlus::StateStruct Gen1State;
		Gen1State.ConvertVecttoState(IArgs.GenState);
		Gen1.setstate(Gen1State);

		Irand.configure(Gen1, XorShiftPlus(), alpha);	// second generator is dummy.
		if (IArgs.Irand.istrulyempty())
			Irand.resize(N);
		else if (IArgs.Irand.size() == N)
			Irand.assign(IArgs.Irand);				// initializing Vector
		else{
			// GIVE ERROR MESSAGE HERE
			return;
		}
		
		int LoopLimit = (8192 <= onemsbyTstep*NoOfms) ? 8192 : onemsbyTstep*NoOfms + 1;
		RandMat[0] = Irand;
		Irand.generator1().getstate().ConvertStatetoVect(GenMat[0]);

		for (int j = 1; j < LoopLimit; ++j){
			Irand.generate();
			RandMat[j] = Irand;
			Irand.generator1().getstate().ConvertStatetoVect(GenMat[j]);
		}
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
	}
	void DoOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars){
		DoFullOutput(StateOut, OutVars);
		if (StorageStepSize && !(Time % (StorageStepSize*onemsbyTstep))){
			DoSparseOutput(StateOut, OutVars);
		}
	}
	void DoSingleStateOutput(SingleStateStruct &FinalStateOut);
private:
	void DoSparseOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars);
	void DoFullOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars);
};

struct OutputVarsStruct{
	MexMatrix<float> WeightOut;
	MexMatrix<float> Iin;
	MexMatrix<float> Itot;
	struct SpikeListStruct{
		MexVector<int> SpikeSynInds;
		MexVector<int> TimeRchdStartInds;
		SpikeListStruct() : SpikeSynInds(), TimeRchdStartInds(){}
	} SpikeList;

	OutputVarsStruct() :
		WeightOut(),
		Itot(),
		Iin(),
		SpikeList(){}

	void initialize(const InternalVars &);
};

struct StateVarsOutStruct{
	MexMatrix<float> WeightOut;
	MexMatrix<float> VOut;
	MexMatrix<float> UOut;
	MexMatrix<float> Iin1Out;
	MexMatrix<float> Iin2Out;
	MexMatrix<float> WeightDerivOut;
	MexMatrix<uint32_t> GenStateOut;
	MexMatrix<float> IrandOut;

	MexVector<int> TimeOut;
	MexVector<MexVector<MexVector<int> > > SpikeQueueOut;
	MexVector<int> CurrentQIndexOut;
	MexMatrix<int> LSTNeuronOut;
	MexMatrix<int> LSTSynOut;

	StateVarsOutStruct() :
		WeightOut(),
		VOut(),
		UOut(),
		Iin1Out(),
		Iin2Out(),
		WeightDerivOut(),
		GenStateOut(),
		IrandOut(),
		TimeOut(),
		SpikeQueueOut(),
		CurrentQIndexOut(),
		LSTNeuronOut(),
		LSTSynOut() {}

	void initialize(const InternalVars &);
};
struct SingleStateStruct{
	MexVector<float> Weight;
	MexVector<float> V;
	MexVector<float> U;
	MexVector<float> Iin1;
	MexVector<float> Iin2;
	MexVector<float> WeightDeriv;
	MexVector<uint32_t> GenState;
	MexVector<float> Irand;

	MexVector<MexVector<int > > SpikeQueue;
	MexVector<int> LSTNeuron;
	MexVector<int> LSTSyn;
	int Time;
	int CurrentQIndex;

	SingleStateStruct() :
		Weight(),
		V(),
		U(),
		Iin1(),
		Iin2(),
		WeightDeriv(),
		GenState(),
		Irand(),
		SpikeQueue(),
		LSTNeuron(),
		LSTSyn() {}

	virtual void initialize(const InternalVars &) {}
};
struct FinalStateStruct : public SingleStateStruct{
	void initialize(const InternalVars &);
};

struct InitialStateStruct : public SingleStateStruct{
	void initialize(const InternalVars &);
};

void CountingSort(int N, MexVector<Synapse> &Network, MexVector<int> &indirection);

void SimulateParallel(
	InputArgs &&InputArguments,
	OutputVarsStruct &PureOutputs,
	StateVarsOutStruct &StateVarsOutput,
	FinalStateStruct &FinalStateOutput,
	InitialStateStruct &InitalStateOutput);

#endif