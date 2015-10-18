#include <vector>
#include <iostream>
#include <tbb\parallel_for.h>
#include <tbb\blocked_range.h>
#include <tbb\atomic.h>
#include <fstream>
#include <chrono>
#include <cmath>

#if defined TIME_DEL_NET_SIM_AS_SUB
	#define HEADER_PATHS_TDNS ..
#elif !defined HEADER_PATHS_TDNS
	#define HEADER_PATHS_TDNS .
#endif

#define SETQUOTE(A) #A
#define JOIN_STRING(A,B,C) SETQUOTE(A##B##C)
#define JOIN_LIB_PATH(PRE, CENT, POST) JOIN_STRING(PRE, CENT, POST)

#include "..\Headers\Network.hpp"
#include "..\Headers\NeuronSim.hpp"

#include "..\Headers\IExtHeaders\IExtCode.hpp"

#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_TDNS, \MexMemoryInterfacing\Headers\MexMem.hpp)
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_TDNS, \MexMemoryInterfacing\Headers\GenericMexIO.hpp)
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_TDNS, \MexMemoryInterfacing\Headers\InterruptHandling.hpp)
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_TDNS, \MexMemoryInterfacing\Headers\LambdaToFunction.hpp)
#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_TDNS, \MexMemoryInterfacing\Headers\FlatVectTree\FlatVectTree.hpp)

#include JOIN_LIB_PATH(..\..\, HEADER_PATHS_TDNS, \RandomNumGen\Headers\FiltRandomTBB.hpp)

#include <emmintrin.h>
#include <smmintrin.h>

using namespace std;

#define VECTOR_IMPLEMENTATION

void CountingSort(int N, MexVector<Synapse> &Network, MexVector<size_t> &indirection)
{
	MexVector<int> CumulativeCountStart(N,0);
	MexVector<int> IOutsertIndex(N);
	size_t M = Network.size();
	if (indirection.size() != M) indirection.resize(M);

	for (int i = 0; i < M; ++i){
		CumulativeCountStart[Network[i].NEnd - 1]++;
	}
	IOutsertIndex[0] = 0;
	for (int i = 1; i < N; ++i){
		CumulativeCountStart[i] += CumulativeCountStart[i - 1];
		IOutsertIndex[i] = CumulativeCountStart[i - 1];
	}
	int indirIter = 0;
	for (int i = 0; i < M; ++i){
		int NeuronEnd = Network[i].NEnd;
		indirection[IOutsertIndex[NeuronEnd-1]] = i;
		IOutsertIndex[NeuronEnd - 1]++;
	}
}

void CurrentUpdate::operator () (const tbb::blocked_range<int*> &BlockedRange) const{
	const float &I0 = IntVars.I0;
	auto &Network           = IntVars.Network;
	auto &Iin1              = IntVars.Iin1;
	auto &Iin2              = IntVars.Iin2;
	auto &LastSpikedTimeSyn = IntVars.LSTSyn;
	auto &time              = IntVars.Time;

	int *begin = BlockedRange.begin();
	int *end = BlockedRange.end();
	for (int * iter = begin; iter < end; ++iter){
		Synapse CurrentSynapse = Network[*iter];
		int CurrentSynapseInd = *iter;
		__m128 AddedCurrent;
		AddedCurrent.m128_u64[0] = 0; AddedCurrent.m128_u64[1] = 0;
		AddedCurrent.m128_f32[0] = CurrentSynapse.Weight;
		AddedCurrent.m128_u32[0] += (32 << 23);

		Iin1[CurrentSynapse.NEnd - 1].fetch_and_add((long long)AddedCurrent.m128_f32[0]);
		Iin2[CurrentSynapse.NEnd - 1].fetch_and_add((long long)AddedCurrent.m128_f32[0]);
		LastSpikedTimeSyn[CurrentSynapseInd] = time;
	}
}

void NeuronSimulate::operator() (tbb::blocked_range<int> &Range) const{

	auto &Vnow = IntVars.V;
	auto &Unow = IntVars.U;

	auto &Neurons = IntVars.Neurons;
	auto &Network = IntVars.Network;

	auto &Iin1 = IntVars.Iin1;
	auto &Iin2 = IntVars.Iin2;
	auto &Iext = IntVars.IextInterface.Iext;

	auto &PreSynNeuronSectionBeg = IntVars.PreSynNeuronSectionBeg;
	auto &PreSynNeuronSectionEnd = IntVars.PreSynNeuronSectionEnd;

	auto &LastSpikedTimeNeuron = IntVars.LSTNeuron;
	auto &onemsbyTstep = IntVars.onemsbyTstep;
	auto &time = IntVars.Time;

	size_t QueueSize = onemsbyTstep*IntVars.DelayRange;
	size_t RangeBeg = Range.begin();
	size_t RangeEnd = Range.end();
	for (size_t j = RangeBeg; j < RangeEnd; ++j){
		if (Vnow[j] == 30.0f){
			//Implementing Izhikevich resetting
			Vnow[j] = Neurons[j].c;
			Unow[j] += Neurons[j].d;
		}
		else{
			//Implementing Izhikevich differential equation
			float Vnew, Unew;
			Vnew = Vnow[j] + (Vnow[j] * (0.04f*Vnow[j] + 5.0f) + 140.0f - Unow[j] + (float)(Iin2[j] - Iin1[j]) / (1i64 << 32) + Iext[j]) / onemsbyTstep;
			Unew = Unow[j] + (Neurons[j].a*(Neurons[j].b*Vnow[j] - Unow[j])) / onemsbyTstep;
			Vnow[j] = (Vnew > -100)? Vnew: -100;
			Unow[j] = Unew;

			//Implementing Network Computation in case a Neuron has spiked in the current interval
			if (Vnow[j] >= 30.0f){
				Vnow[j] = 30.0f;
				//NSpikesGenminProc += ((PreSynNeuronSectionBeg[j] >= 0) ? PreSynNeuronSectionEnd[j] - PreSynNeuronSectionBeg[j] : 0);
				LastSpikedTimeNeuron[j] = time;
				//Space to implement any causal Learning Rule
			}
		}
	}
}
void CurrentAttenuate::operator() (tbb::blocked_range<int> &Range) const {

	auto &Iin1 = IntVars.Iin1;
	auto &Iin2 = IntVars.Iin2;
	auto &attenFactor1 = IntVars.CurrentDecayFactor1;
	auto &attenFactor2 = IntVars.CurrentDecayFactor2;

	tbb::atomic<long long> *Begin1 = &Iin1[Range.begin()];
	tbb::atomic<long long> *End1 = &Iin1[Range.end()-1] + 1;
	tbb::atomic<long long> *Begin2 = &Iin2[Range.begin()];
	tbb::atomic<long long> *End2 = &Iin2[Range.end() - 1] + 1;

	for (tbb::atomic<long long> *i = Begin1, *j = Begin2; i < End1; ++i, ++j){
		(*i) = (long long)(float(i->load()) * attenFactor1);
		(*j) = (long long)(float(j->load()) * attenFactor2);
	}
}

void StateVarsOutStruct::initialize(const InternalVars &IntVars) {

	auto onemsbyTstep = IntVars.onemsbyTstep;
	auto NoOfms = IntVars.NoOfms;
	auto StorageStepSize = IntVars.StorageStepSize;
	auto Tbeg = IntVars.Time;
	auto nSteps = IntVars.nSteps;
	auto OutputControl = IntVars.OutputControl;
	auto N = IntVars.N;
	auto M = IntVars.M;
	auto DelayRange = IntVars.DelayRange;

	if (OutputControl & OutOps::WEIGHT_REQ)
		if (!(IntVars.InterestingSyns.size()))
			this->WeightOut = MexMatrix<float>(0, M);

	if (OutputControl & OutOps::V_REQ)
		this->VOut = MexMatrix<float>(0, N);

	if (OutputControl & OutOps::U_REQ)
		this->UOut = MexMatrix<float>(0, N);

	if (OutputControl & OutOps::I_IN_1_REQ)
		this->Iin1Out = MexMatrix<float>(0, N);
	
	if (OutputControl & OutOps::I_IN_2_REQ)
		this->Iin2Out = MexMatrix<float>(0, N);

	if (OutputControl & OutOps::WEIGHT_DERIV_REQ)
		this->WeightDerivOut = MexMatrix<float>(0, M);

	this->IextInterface.initialize(IntVars.IextInterface, IntVars);

	this->TimeOut = MexVector<int>(0);

	if (OutputControl & OutOps::LASTSPIKED_NEU_REQ)
		this->LSTNeuronOut = MexMatrix<int>(0, N);

	if (OutputControl & OutOps::LASTSPIKED_SYN_REQ)
		this->LSTSynOut = MexMatrix<int>(0, M);

	if (OutputControl & OutOps::SPIKE_QUEUE_REQ)
		this->SpikeQueueOut = FlatVectTree<int>(2);

	if (OutputControl & OutOps::CURRENT_QINDS_REQ)
		this->CurrentQIndexOut = MexVector<int>(0);
}
void OutputVarsStruct::initialize(const InternalVars &IntVars){
	size_t TimeDimLen;
	auto N = IntVars.N;
	auto onemsbyTstep = IntVars.onemsbyTstep;
	auto NoOfms = IntVars.NoOfms;
	auto StorageStepSize = IntVars.StorageStepSize;
	auto Tbeg = IntVars.Time;
	auto nSteps = IntVars.nSteps;
	auto OutputControl = IntVars.OutputControl;

	if (OutputControl & OutOps::WEIGHT_REQ)
		if (IntVars.InterestingSyns.size())
			this->WeightOut = MexMatrix<float>(0, IntVars.InterestingSyns.size());
	if (OutputControl & OutOps::I_IN_REQ)
		this->Iin = MexMatrix<float>(0, N);
	if (OutputControl & OutOps::I_TOT_REQ)
		this->Itot = MexMatrix<float>(0, N);
	this->NoOfSpikes = MexVector<uint64_t>(0);

	// Initializing Output Variables for IextInterface
	this->IextInterface.initialize(IntVars.IextInterface, IntVars);

	if ((OutputControl & OutOps::SPIKE_LIST_REQ) && !StorageStepSize);
		// The vector is initialized to size zero regardless. the If condition is 
		// just kept for code conformity
}
void SingleStateStruct::initialize(const InternalVars &IntVars){
	auto OutputControl	= IntVars.OutputControl;
	auto DelayRange		= IntVars.DelayRange;
	auto onemsbyTstep	= IntVars.onemsbyTstep;
	auto N				= IntVars.N;
	auto M				= IntVars.M;

	this->V = MexVector<float>(N);
	this->U = MexVector<float>(N);
	this->Iin1 = MexVector<float>(N);
	this->Iin2 = MexVector<float>(N);
	this->WeightDeriv = MexVector<float>(M);
	this->IextInterface.initialize(IntVars.IextInterface, IntVars);
	this->Weight = MexVector<float>(M);
	this->LSTNeuron = MexVector<int>(N);
	this->LSTSyn = MexVector<int>(M);
	this->SpikeQueue = FlatVectTree<int>(1);
	this->CurrentQIndex = -1;
	this->Time = -1;
}
void InternalVars::DoOutput(StateVarsOutStruct &StateOut, OutputVarsStruct &OutVars) {

	if (StorageStepSize && (Time % (StorageStepSize*onemsbyTstep) == 0) || !StorageStepSize) {
		
		// Storing U,V,Iin and Time
		if (OutputControl & OutOps::V_REQ)
			StateOut.VOut.push_row(V);
		if (OutputControl & OutOps::U_REQ)
			StateOut.UOut.push_row(U);
		if (OutputControl & OutOps::I_IN_1_REQ) {
			StateOut.Iin1Out.push_row_size(1);
			for (int j = 0; j < N; ++j)
				StateOut.Iin1Out.lastRow()[j] = (float)Iin1[j] / (1i64 << 32);
		}
		if (OutputControl & OutOps::I_IN_2_REQ) {
			StateOut.Iin2Out.push_row_size(1);
			for (int j = 0; j < N; ++j)
				StateOut.Iin2Out.lastRow()[j] = (float)Iin2[j] / (1i64 << 32);
		}

		//Storing Weight Derivative
		if (OutputControl & OutOps::WEIGHT_DERIV_REQ) {
			StateOut.WeightDerivOut.push_row(WeightDeriv);
		}

		// Storing IextInterface related state / output vars
		IExtInterface::doOutput(
			StateOut.IextInterface,
			OutVars.IextInterface,
			this->IextInterface,
			*this
			);

		// Storing Current Time Instant
		StateOut.TimeOut.push_back(Time);

		// Storing Weights
		if (OutputControl & OutOps::WEIGHT_REQ && InterestingSyns.size()) {
			OutVars.WeightOut.push_row_size(1);
			for (int j = 0; j < InterestingSyns.size(); ++j)
				OutVars.WeightOut.lastRow()[j] = Network[InterestingSyns[j]].Weight;
		}
		else if (OutputControl & OutOps::WEIGHT_REQ) {
			StateOut.WeightOut.push_row_size(1);
			for (int j = 0; j < M; ++j)
				StateOut.WeightOut.lastRow()[j] = Network[j].Weight;
		}

		// Storing Spike Queue related state informations
		if (OutputControl & OutOps::SPIKE_QUEUE_REQ)
			StateOut.SpikeQueueOut.push_back(SpikeQueue);
		if (OutputControl & OutOps::CURRENT_QINDS_REQ)
			StateOut.CurrentQIndexOut.push_back(CurrentQIndex);

		// Storing last Spiked timings
		if (OutputControl & OutOps::LASTSPIKED_NEU_REQ)
			StateOut.LSTNeuronOut.push_row(LSTNeuron);
		if (OutputControl & OutOps::LASTSPIKED_SYN_REQ)
			StateOut.LSTSynOut.push_row(LSTSyn);

		// Storing Iin
		if (OutputControl & OutOps::I_IN_REQ) {
			OutVars.Iin.push_row_size(1);
			for (int j = 0; j < N; ++j)
				OutVars.Iin.lastRow()[j] = (float)(Iin2[j] - Iin1[j]) / (1i64 << 32);
		}

		// Storing Itot
		if (OutputControl & OutOps::I_TOT_REQ) {
			OutVars.Itot.push_row_size(1);
			for (int j = 0; j < N; ++j)
				OutVars.Itot.lastRow()[j] = this->IextInterface.Iext[j] + (float)(Iin2[j] - Iin1[j]) / (1i64 << 32);
		}

		// Storing NoOfSpikes
		OutVars.NoOfSpikes.push_back(NoOfSpikes);
		NoOfSpikes = 0;
	}

	// Storing Spike List (Only if StorageStepSize == 0
	if (OutputControl & OutOps::SPIKE_LIST_REQ) {
		OutVars.SpikeList.TimeRchd         .push_back(Time+1);
		OutVars.SpikeList.TimeRchdStartInds.push_back(OutVars.SpikeList.SpikeSynInds.size());
		for (auto Spike : SpikeQueue[CurrentQIndex]) {
			OutVars.SpikeList.SpikeSynInds.push_back(Spike);
		}
		if (i == nSteps) {
			// Storing spikes which are generated but not gonna arrive next turn
			for (int j = 1; j < DelayRange*onemsbyTstep; ++j) {
				OutVars.SpikeList.TimeRchd         .push_back(Time + j + 1);
				OutVars.SpikeList.TimeRchdStartInds.push_back(OutVars.SpikeList.SpikeSynInds.size());
				for (auto Spike : SpikeQueue[(CurrentQIndex + j) % (onemsbyTstep*DelayRange)]) {
					OutVars.SpikeList.SpikeSynInds.push_back(Spike);
				}
			}
			OutVars.SpikeList.TimeRchdStartInds.push_back(OutVars.SpikeList.SpikeSynInds.size());
			// Final push_back in order to be able to infer end.
		}
	}
}
void InternalVars::DoSingleStateOutput(SingleStateStruct &SingleStateOut){
	size_t QueueSize = onemsbyTstep * DelayRange;
	for (int j = 0; j < N; ++j){
		SingleStateOut.Iin1[j] = (float)Iin1[j] / (1i64 << 32);
		SingleStateOut.Iin2[j] = (float)Iin2[j] / (1i64 << 32);
	}
	// Storing IextInterface related state vars
	IExtInterface::doSingleStateOutput(
		SingleStateOut.IextInterface,
		this->IextInterface,
		*this
		);

	SingleStateOut.V = V;
	SingleStateOut.U = U;
	for (int j = 0; j < M; ++j){
		SingleStateOut.Weight[j] = Network[j].Weight;
	}
	SingleStateOut.WeightDeriv = WeightDeriv;
	SingleStateOut.SpikeQueue.append(SpikeQueue);
	SingleStateOut.CurrentQIndex = CurrentQIndex;
	SingleStateOut.LSTNeuron = LSTNeuron;
	SingleStateOut.LSTSyn = LSTSyn;
	SingleStateOut.Time = Time;
}
void InternalVars::DoInputStateOutput(InputArgs &InputStateOut){
	
	// Input Vectors
	InputStateOut.NStart.resize(M);
	InputStateOut.NEnd  .resize(M);
	InputStateOut.Weight.resize(M);
	InputStateOut.Delay .resize(M);
	
	MexTransform(Network.begin(), Network.end(), InputStateOut.NStart.begin(), FFL([ ](Synapse &Syn)->int  {return Syn.NStart       ; }));
	MexTransform(Network.begin(), Network.end(), InputStateOut.NEnd  .begin(), FFL([ ](Synapse &Syn)->int  {return Syn.NEnd         ; }));
	MexTransform(Network.begin(), Network.end(), InputStateOut.Weight.begin(), FFL([ ](Synapse &Syn)->float{return Syn.Weight       ; }));
	MexTransform(Network.begin(), Network.end(), InputStateOut.Delay .begin(), FFL([&](Synapse &Syn)->float{return (float)Syn.DelayinTsteps / onemsbyTstep; }));

	InputStateOut.a.resize(N);
	InputStateOut.b.resize(N);
	InputStateOut.c.resize(N);
	InputStateOut.d.resize(N);

	MexTransform(Neurons.begin(), Neurons.end(), InputStateOut.a.begin(), FFL([](Neuron &Neu)->float{return Neu.a; }));
	MexTransform(Neurons.begin(), Neurons.end(), InputStateOut.b.begin(), FFL([](Neuron &Neu)->float{return Neu.b; }));
	MexTransform(Neurons.begin(), Neurons.end(), InputStateOut.c.begin(), FFL([](Neuron &Neu)->float{return Neu.c; }));
	MexTransform(Neurons.begin(), Neurons.end(), InputStateOut.d.begin(), FFL([](Neuron &Neu)->float{return Neu.d; }));

	InputStateOut.InterestingSyns = InterestingSyns;

	// Compulsory Simulation Parameters
	InputStateOut.onemsbyTstep = onemsbyTstep;
	InputStateOut.NoOfms       = NoOfms;
	InputStateOut.DelayRange   = DelayRange;

	// Optional Simulation Parameters
	InputStateOut.OutputControlString   = OutputControlString   ;
	InputStateOut.OutputControl         = OutputControl         ;
	InputStateOut.StorageStepSize       = StorageStepSize       ;
	InputStateOut.StatusDisplayInterval = StatusDisplayInterval ;

	// Optional Simulation Algorithm Parameters
	InputStateOut.I0                  = I0                  ;
	InputStateOut.CurrentDecayFactor1 = CurrentDecayFactor1 ;
	InputStateOut.CurrentDecayFactor2 = CurrentDecayFactor2 ;

	// Assigining IExt Input Vars
	IExtInterface::doInputVarsOutput(
		InputStateOut.IextInterface,
		this->IextInterface,
		*this);

	// Initializing and Assigning the State Variables
	InputStateOut.InitialState.initialize(*this);
	DoSingleStateOutput(InputStateOut.InitialState);
	IExtInterface::doSingleStateOutput(
		InputStateOut.InitialState.IextInterface, 
		this->IextInterface, 
		*this);
}

void CachedSpikeStorage(InternalVars &IntVars){

	auto &N = IntVars.N;
	auto &Network = IntVars.Network;

	auto &SpikeQueue = IntVars.SpikeQueue;
	auto &LastSpikedTimeNeuron = IntVars.LSTNeuron;

	auto &preSynNeuronSectionBeg = IntVars.PreSynNeuronSectionBeg;
	auto &preSynNeuronSectionEnd = IntVars.PreSynNeuronSectionEnd;

	auto &CurrentQIndex = IntVars.CurrentQIndex;
	auto &time = IntVars.Time;

	const size_t nBins = IntVars.onemsbyTstep * IntVars.DelayRange;
	const size_t CacheBuffering = 128;	// Each time a cache of size 64 will be pulled in 
	
	auto &BinningBuffer = IntVars.BinningBuffer;
	auto &BufferInsertIndex = IntVars.BufferInsertIndex;
	auto &AddressOffset = IntVars.AddressOffset;
	
	for (int j = 0; j < nBins; ++j){
		AddressOffset[j] = (reinterpret_cast<size_t>(SpikeQueue[j].end()) & 0x0F) >> 2;
	}
	for (int j = 0; j < N; ++j){
		if (LastSpikedTimeNeuron[j] == time){
			size_t k = preSynNeuronSectionBeg[j];
			size_t kend = preSynNeuronSectionEnd[j];

			if (k != kend){
				size_t NoofCurrNeuronSpikes = kend - k;
				MexVector<Synapse>::iterator iSyn = Network.begin() + k;
				MexVector<Synapse>::iterator iSynEnd = Network.begin() + kend;

				//TotalSpikesTemp += iSynEnd - iSyn;
				for (; iSyn < iSynEnd; ++iSyn, ++k){
					int CurrIndex = (CurrentQIndex + iSyn->DelayinTsteps) % nBins;
					int CurrAddressOffset = AddressOffset[CurrIndex];
					int BufferIndex = BufferInsertIndex[CurrIndex];
					int *BufferIndexPtr = &BufferInsertIndex[CurrIndex];
					int *CurrAddressOffsetPtr = &AddressOffset[CurrIndex];

					if (BufferIndex == CacheBuffering){
						if (CurrAddressOffset){
							for (int k = 0; k < 4-CurrAddressOffset; ++k){
								SpikeQueue[CurrIndex].push_back(*(reinterpret_cast<int*>(BinningBuffer.begin() + (CurrIndex + 1)*CacheBuffering / 4)
									                                      - (4-CurrAddressOffset) + k));
								// This is bcuz, the buffer at CurrIndex ends at
								// BinningBuffer.begin() + (CurrIndex + 1)*CacheBuffering / 4 [__m128*]
								// therefore we move 4-CurrAddressOffset (which is the no. of elems to be inserted)
								// behind and push back.
								}
							BufferIndex -= 4 - CurrAddressOffset;
							*BufferIndexPtr -= 4 - CurrAddressOffset;
							*CurrAddressOffsetPtr = 0;
						}
						else{
							SpikeQueue[CurrIndex].push_size(CacheBuffering);
							//TotalSpikesTemp += CacheBuffering;
							__m128i* kbeg = reinterpret_cast<__m128i*>(SpikeQueue[CurrIndex].end() - CacheBuffering);
							__m128i* kend = reinterpret_cast<__m128i*>(SpikeQueue[CurrIndex].end());
							__m128i* lbeg = reinterpret_cast<__m128i*>(BinningBuffer.begin() + CurrIndex * CacheBuffering / 4);

							for (__m128i* k = kbeg, *l = lbeg; k < kend; k++, l++){
								_mm_stream_si128(k, _mm_load_si128(l));
							}
							BufferIndex = 0;
							*BufferIndexPtr = 0;
						}
						
					}

					reinterpret_cast<int *>(BinningBuffer.begin())[CurrIndex*CacheBuffering + BufferIndex]
						= k;
					++*BufferIndexPtr;
				}
			}
		}
	}

	for (int i = 0; i < nBins; ++i){
		size_t CurrNElems = BufferInsertIndex[i];
		SpikeQueue[i].push_size(CurrNElems);
		int* kbeg = SpikeQueue[i].end() - CurrNElems;
		int* kend = SpikeQueue[i].end();
		int* lbeg = reinterpret_cast<int*>(BinningBuffer.begin() + i * CacheBuffering/4);
		for (int* k = kbeg, *l = lbeg; k < kend; ++k, ++l){
			*k = *l;
		}
		BufferInsertIndex[i] = 0;
	}
}

void SimulateParallel(
	InputArgs &&InputArguments,
	OutputVarsStruct &PureOutputs,
	StateVarsOutStruct &StateVarsOutput,
	SingleStateStruct &FinalStateOutput,
	InputArgs &InputStateOutput
)
{
	// Aliasing Input Arguments Into Appropriate
	// "In Function" state and input variables

	// Initialization and aliasing of All the input / State / parameter variables.
	InternalVars IntVars(InputArguments);

	// Aliasing of Data members in IntVar
	MexVector<Synapse>			&Network              = IntVars.Network;
	MexVector<Neuron>			&Neurons              = IntVars.Neurons;
	MexVector<float>			&Vnow                 = IntVars.V;
	MexVector<float>			&Unow                 = IntVars.U;
	MexVector<int>				&InterestingSyns      = IntVars.InterestingSyns;
	atomicLongVect				&Iin1                 = IntVars.Iin1;
	atomicLongVect				&Iin2                 = IntVars.Iin2;
	MexVector<MexVector<int> >	&SpikeQueue           = IntVars.SpikeQueue;
	MexVector<int>				&LastSpikedTimeNeuron = IntVars.LSTNeuron;
	MexVector<int>				&LastSpikedTimeSyn    = IntVars.LSTSyn;
	

	size_t &NoOfms              = IntVars.NoOfms;
	const size_t &onemsbyTstep  = IntVars.onemsbyTstep;
	size_t Tbeg                 = IntVars.Time;				//Tbeg is just an initial constant, 
	size_t &time                = IntVars.Time;				//time is the actual changing state variable
	size_t &DelayRange          = IntVars.DelayRange;
	size_t &CurrentQueueIndex   = IntVars.CurrentQIndex;
	size_t &StorageStepSize     = IntVars.StorageStepSize;
	size_t &OutputControl       = IntVars.OutputControl;
	size_t &i                   = IntVars.i;
	size_t &nSteps              = IntVars.nSteps;
	uint64_t &NoOfSpikes        = IntVars.NoOfSpikes;

	const float &I0			= IntVars.I0;	// Value of the current factor to be multd with weights (constant)
	// calculate value of alpha for filtering
	// alpha = 0 => no filtering
	// alpha = 1 => complete filtering
	const float &CurrentDecayFactor1	= IntVars.CurrentDecayFactor1;	//Current Decay Factor in the current model (possibly input in future)
	const float &CurrentDecayFactor2	= IntVars.CurrentDecayFactor2;
	const size_t &StatusDisplayInterval = IntVars.StatusDisplayInterval;

	// other data members. probably derived from inputs or something
	// I think should be a constant. (note that it is possible that 
	// I club some of these with the inputs in future revisions like
	// CurrentDecayFactor

	size_t QueueSize = SpikeQueue.size();
	size_t N = IntVars.Neurons.size(), M = IntVars.Network.size();			
	
	// VARIOuS ARRAYS USED apart from those in the argument list and Output List.
	// Id like to call them intermediate arrays, required for simulation but are
	// not state, input or output vectors.
	// they are typically of the form of some processed version of an input vector
	// thus they dont change with time and are prima facie not used to generate output

	MexVector<size_t> &AuxArray                = IntVars.AuxArray;

	MexVector<size_t> &PreSynNeuronSectionBeg  = IntVars.PreSynNeuronSectionBeg;	
	
	MexVector<size_t> &PreSynNeuronSectionEnd  = IntVars.PreSynNeuronSectionEnd;

	MexVector<size_t> &PostSynNeuronSectionBeg = IntVars.PostSynNeuronSectionBeg;
	
	MexVector<size_t> &PostSynNeuronSectionEnd = IntVars.PostSynNeuronSectionEnd;
	
	
	//----------------------------------------------------------------------------------------------//
	//--------------------------------- Initializing output Arrays ---------------------------------//
	//----------------------------------------------------------------------------------------------//

	StateVarsOutput.initialize(IntVars);
	PureOutputs.initialize(IntVars);
	if (OutputControl & OutOps::FINAL_STATE_REQ) FinalStateOutput.initialize(IntVars);
	
	//---------------------------- Initializing the Intermediate Arrays ----------------------------//
	CountingSort(N, Network, AuxArray);	// Perform counting sort by (NEnd, NStart)
	                                    // to get AuxArray

	// Sectioning the Network and AuxArray Arrays as according to 
	// definition of respective variables above
	PreSynNeuronSectionBeg[Network[0].NStart - 1] = 0;
	PostSynNeuronSectionBeg[Network[AuxArray[0]].NEnd - 1] = 0;

	PreSynNeuronSectionBeg[0] = 0;
	PreSynNeuronSectionEnd[0] = 0;
	PostSynNeuronSectionBeg[0] = 0;
	PostSynNeuronSectionEnd[0] = 0;

	PreSynNeuronSectionEnd[Network[M - 1].NStart - 1] = M;
	PostSynNeuronSectionEnd[Network[AuxArray[M - 1]].NEnd - 1] = M;

	for (int j = 1; j<M; ++j){
		if (Network[j - 1].NStart != Network[j].NStart){
			PreSynNeuronSectionBeg[Network[j].NStart - 1] = j;
			PreSynNeuronSectionEnd[Network[j - 1].NStart - 1] = j;
		}
		if (Network[AuxArray[j - 1]].NEnd != Network[AuxArray[j]].NEnd){
			PostSynNeuronSectionBeg[Network[AuxArray[j]].NEnd - 1] = j;
			PostSynNeuronSectionEnd[Network[AuxArray[j - 1]].NEnd - 1] = j;
		}
	}
	for (int j = 1; j < N; ++j){
		if (PreSynNeuronSectionBeg[j] == -1){
			PreSynNeuronSectionBeg[j] = PreSynNeuronSectionEnd[j]
				= PreSynNeuronSectionEnd[j - 1];
		}
		if (PostSynNeuronSectionBeg[j] == -1){
			PostSynNeuronSectionBeg[j] = PostSynNeuronSectionEnd[j]
				= PostSynNeuronSectionEnd[j - 1];
		}
	}

	// Here I assume that the first neuron is always excitatory
	// and all excitatory neurons are stored contiguously starting 
	// from the first neuron, and that excitatory and inhibitory
	// neurons differe in their parameter a (which is tru in izhikevich case)

	// The Structure of iteration below is given below
	/*
		->Iin[t-1] ----------
		|                    \
		|                     \
		|                      >Itemp ---------- * CurrentDecayFactor ---- Iin[t]
		|                     /             \                                |
		|                    /               \                               |
		->SpikeQueue[t-1] ---                 > V,U[t] ----- SpikeQueue[t]   |
		|       |                            /    |                |         |
		|       |               V,U[t-1] ----     |                |         |
		|       |                  |              |                |         |
		|===<===|====<=========<===|==Loop Iter<==|=======<========|===<=====|

		The vector corresponding to the spikes processed in the current 
		iteration is cleared after the calculation of Itemp
	*/
	// Giving Initial State if Asked For
	
	if (OutputControl & OutOps::INITIAL_STATE_REQ){
		IntVars.DoInputStateOutput(InputStateOutput);
	}
	size_t maxSpikeno = 0;
	tbb::affinity_partitioner apCurrentUpdate;
	tbb::affinity_partitioner apNeuronSim;
	size_t TotalStorageStepSize = (StorageStepSize*onemsbyTstep); // used everywhere
	int epilepsyctr = 0;

	// ------------------------------------------------------------------------------ //
	// ------------------------------ Simulation Loop ------------------------------- //
	// ------------------------------------------------------------------------------ //
	// Profiling Times.
	size_t IExtGenTime = 0,
		   IRandGenTime = 0,
		   IUpdateTime = 0,
		   SpikeStoreTime = 0,
		   NeuronCalcTime = 0,
		   OutputTime = 0;
	std::chrono::system_clock::time_point
		IExtGenTimeBeg, IExtGenTimeEnd,
		IRandGenTimeBeg, IRandGenTimeEnd,
		IUpdateTimeBeg, IUpdateTimeEnd,
		SpikeStoreTimeBeg, SpikeStoreTimeEnd,
		NeuronCalcTimeBeg, NeuronCalcTimeEnd,
		OutputTimeBeg, OutputTimeEnd;

	for (i = 1; i<=nSteps; ++i){
		
		time = time + 1;

		// Updating Iext for current time instant
		IExtGenTimeBeg = std::chrono::system_clock::now();
		IExtInterface::updateIExt(IntVars.IextInterface, IntVars);
		IExtGenTimeEnd = std::chrono::system_clock::now();
		IExtGenTime += std::chrono::duration_cast<std::chrono::microseconds>(IExtGenTimeEnd - IExtGenTimeBeg).count();

		// This iteration applies time update equation for internal current
		// in this case, it is just an exponential attenuation
		tbb::parallel_for(tbb::blocked_range<int>(0, N, 3000),
			CurrentAttenuate(IntVars));

		size_t QueueSubEnd = SpikeQueue[CurrentQueueIndex].size();

		IUpdateTimeBeg = std::chrono::system_clock::now();
		// This iter calculates Itemp as in above diagram
		if (SpikeQueue[CurrentQueueIndex].size() != 0)
			tbb::parallel_for(tbb::blocked_range<int*>((int*)&SpikeQueue[CurrentQueueIndex][0],
				(int*)&SpikeQueue[CurrentQueueIndex][QueueSubEnd - 1] + 1, 10000), 
				CurrentUpdate(IntVars), apCurrentUpdate);
		SpikeQueue[CurrentQueueIndex].clear();
		IUpdateTimeEnd = std::chrono::system_clock::now();
		IUpdateTime += std::chrono::duration_cast<std::chrono::microseconds>(IUpdateTimeEnd - IUpdateTimeBeg).count();

		// Calculation of V,U[t] from V,U[t-1], Iin = Itemp
		NeuronCalcTimeBeg = std::chrono::system_clock::now();
		tbb::parallel_for(tbb::blocked_range<int>(0, N, 100), NeuronSimulate(IntVars), apNeuronSim);
		NeuronCalcTimeEnd = std::chrono::system_clock::now();
		NeuronCalcTime += std::chrono::duration_cast<std::chrono::microseconds>(NeuronCalcTimeEnd - NeuronCalcTimeBeg).count();

		// This is code for storing spikes
		SpikeStoreTimeBeg = std::chrono::system_clock::now();
		CachedSpikeStorage(IntVars);
		SpikeStoreTimeEnd = std::chrono::system_clock::now();
		SpikeStoreTime += std::chrono::duration_cast<std::chrono::microseconds>(SpikeStoreTimeEnd - SpikeStoreTimeBeg).count();

		CurrentQueueIndex = (CurrentQueueIndex + 1) % QueueSize;

		maxSpikeno += QueueSubEnd; // Added the number of spikes that arrived in the current time instant
		NoOfSpikes += QueueSubEnd;

		// Epilepsy Check
		if (QueueSubEnd > (2*M) / (5)){
			epilepsyctr++;
			if (epilepsyctr > 100){
				nSteps = i; // Make The current Iteration the last loop iteration
				WriteOutput("Epilepsy Nyuh!! Simulation Terminated at the end of time T = %d steps\n", time);
			}
		}

		// Check for Ctrl-C Interrupt
		if (IsProgramInterrupted()) {
			nSteps = i;   // Make The current Iteration the last loop iteration
			WriteOutput("Simulation Aborted at the end of time T = %d steps\n", time);
			ResetInterrupt();
		}

		OutputTimeBeg = std::chrono::system_clock::now();
		IntVars.DoOutput(StateVarsOutput, PureOutputs);
		OutputTimeEnd = std::chrono::system_clock::now();
		OutputTime += std::chrono::duration_cast<std::chrono::microseconds>(OutputTimeEnd - OutputTimeBeg).count();

		// Status Display Section
		if (!(time % StatusDisplayInterval)){
			WriteOutput("Completed  %d steps with Total no. of Spikes = %d\n", time / (1000 * onemsbyTstep), maxSpikeno);
			maxSpikeno = 0;
		}
	}
	// This needs to be done as the state corresponds to i = nSteps not i = nSteps+1
	i = i - 1;
	if ((OutputControl & OutOps::FINAL_STATE_REQ)){
		IntVars.DoSingleStateOutput(FinalStateOutput);
	}

	WriteOutput("CurrentExt Generation Time = %d millisecs\n", IExtGenTime / 1000);
	WriteOutput("CurrentRand Generation Time = %d millisecs\n", IRandGenTime / 1000);
	WriteOutput("Current Update Time = %d millisecs\n", IUpdateTime / 1000);
	WriteOutput("Spike Storage Time = %d millisecs\n", SpikeStoreTime / 1000);
	WriteOutput("Nuron Calculation Time = %d millisecs\n", NeuronCalcTime / 1000);
	WriteOutput("Output Time = %d millisecs\n", OutputTime / 1000);
}
