#include "getRespSpikes.hpp"
#include <functional>
#include <cmath>
#include <MexMemoryInterfacing/Headers/InterruptHandling.hpp>
#include <MexMemoryInterfacing/Headers/FlatVectTree/FlatVectTree.hpp>

void getRespSpikes::SimulationClass::GenSpikeListCalc(void) {
	// The Algorithm calculates GenSpikeList as follows:
	//
	// # Creating GenSpikeList
	// for each TimeRchd in SpikeList:
	//     get CurrArrivingSyns;
	//     for each syn in CurrArrivingSyns:
	//         GenTime = syn.Delay - TimeRchd
	//         if GenTime >= StartTime && GenTime > GenSpikeList[syn.NStart-1].last():
	//             GenSpikeList[syn.NStart-1].push_back(GenTime)

	for(uint32_t i=0; i < T; ++i) {
		auto CurrentSpikeListBeg = SpikeList.TimeRchdStartInds[i] - StartSpikeListIndex;
		auto CurrentSpikeListEnd = SpikeList.TimeRchdStartInds[i+1] - StartSpikeListIndex;
		auto SpikeArrTime = SpikeList.TimeRchd[i]; // Spike ARRival time
		for(auto j = CurrentSpikeListBeg; j < CurrentSpikeListEnd; ++j) {
			auto &CurrentSpikingSyn = Network[SpikeList.SpikeSynInds[j]];
			auto SpikeGenTime = SpikeArrTime - CurrentSpikingSyn.DelayinTsteps;
			auto &CurrGenSpikeList = GenSpikeList[CurrentSpikingSyn.NStart-1];
			if (SpikeGenTime >= StartTime
			    && (CurrGenSpikeList.isempty()
			        || CurrGenSpikeList.last() < SpikeGenTime)) {
				CurrGenSpikeList.push_back(SpikeGenTime);
			}
		}
	}
}

void getRespSpikes::SimulationClass::SpikeValidityCalc(void) {

	GenSpikeListCalc();

	std::function<bool(const float &, const float &)> isInReset = [&] (const float &UNow, const float &VNow) -> bool {
		float Discriminant = 25.0f-4.0f*0.04f*(140.0f-UNow);
		return Discriminant > 0 && VNow < (-5 - std::sqrt(Discriminant))*12.5f + 5.0f;
	};

	for(uint32_t i=0; i<N; ++i) {
		auto NSpikesGen = GenSpikeList[i].size();
		for(uint32_t j=0; j<NSpikesGen; ++j) {
			auto CurrentSpikeTime = GenSpikeList[i][j];
			auto FirstTimeAfterPrevSpike = (j > 0) ? GenSpikeList[i][j-1]+1 : StartTime;

			DiscreteRange CurrentRange;
			CurrentRange.endPoint = CurrentSpikeTime+1;
			bool foundRangeBegin = false;
			for(uint32_t TIndex=CurrentSpikeTime; TIndex --> FirstTimeAfterPrevSpike;) {
				if (isInReset(U(TIndex-StartTime, i), V(TIndex-StartTime, i))) {
					foundRangeBegin = true;
					CurrentRange.beginPoint = TIndex;
					break;
				}
			}
			if (!foundRangeBegin) {
				CurrentRange.beginPoint = FirstTimeAfterPrevSpike;
			}
			SpikeValidity[i].push_back(CurrentRange);
		}
	}
}

void getRespSpikes::SimulationClass::SpikeTippingCalc(void) {

	std::function<bool(const float &, const float &)> isUnTipped = [&] (const float &UNow, const float &VNow) -> bool {
		float Discriminant = 25.0f-4.0f*0.04f*(140.0f-UNow);
		return Discriminant > 0 && VNow < (-5 + std::sqrt(Discriminant))*12.5f;
	};

	for(uint32_t i=0; i<N; ++i) {
		auto NSpikesGen = GenSpikeList[i].size();
		for(uint32_t j=0; j<NSpikesGen; ++j) {
			auto CurrentSpikeTime = GenSpikeList[i][j];
			auto FirstTimeAfterPrevSpike = (j > 0) ? GenSpikeList[i][j-1]+1 : StartTime;

			DiscreteRange CurrentRange;
			bool foundRangeBegin = false;
			for(uint32_t TIndex=CurrentSpikeTime; TIndex --> FirstTimeAfterPrevSpike;) {
				if (isUnTipped(U(TIndex-StartTime, i), V(TIndex-StartTime, i))) {
					foundRangeBegin = true;
					CurrentRange.beginPoint = TIndex+1;
					CurrentRange.endPoint = CurrentRange.beginPoint+1;
					break;
				}
			}
			if (!foundRangeBegin) {
				CurrentRange.beginPoint = 0;
				CurrentRange.endPoint = 0;
			}
			SpikeTipping[i].push_back(CurrentRange);
		}
	}
}

void getRespSpikes::SimulationClass::initialize(const mxArray *InputmxArray) {
	// Parameters that are initialized are:
	//
	// MexMatrix<float> U;
	// MexMatrix<float> V;
	// MexVector<Synapse> Network;
	//
	// SpikeListStruct SpikeList;

	// Initializing N, M, T
	getInputfromStruct<int32_t>(InputmxArray, "N", this->N, getInputOps(1, "is_required"));
	this->M = FieldInfo<MexVector<int32_t>>::getSize(
		getValidStructField<MexVector<int32_t>>(InputmxArray, "NStart", getInputOps(1, "is_required")));
	this->T = FieldInfo<MexMatrix<float>>::getSize(
		getValidStructField<MexMatrix<float>>(InputmxArray, "U", getInputOps(1, "is_required")), 1);

	// Initializing onemsbyTstep
	getInputfromStruct<int32_t>(InputmxArray, "onemsbyTstep", this->onemsbyTstep, getInputOps(1, "required"));
	
	// Initializing NStart, NEnd, Weight, Delay
	// and Network
	MexVector<int32_t> NStart, NEnd;
	MexVector<float> Weight, Delay;
	getInputfromStruct<int32_t>(InputmxArray, "NStart", NStart, getInputOps(2, "is_required", "required_size", M));
	getInputfromStruct<int32_t>(InputmxArray, "NEnd"  , NEnd  , getInputOps(2, "is_required", "required_size", M));
	getInputfromStruct<float>  (InputmxArray, "Weight", Weight, getInputOps(2, "is_required", "required_size", M));
	getInputfromStruct<float>  (InputmxArray, "Delay" , Delay , getInputOps(2, "is_required", "required_size", M));

	this->Network.resize(M);
	for(int i=0; i<M; ++i) {
		this->Network[i].NStart = NStart[i];
		this->Network[i].NEnd   = NEnd[i];
		this->Network[i].Weight = Weight[i];
		this->Network[i].DelayinTsteps = int32_t(Delay[i]*onemsbyTstep + 0.5);
	}

	// Initializing MExc and NExc
	if (Network[0].Weight < 0) MExc = 0;
	else for(MExc=1; MExc < M && Network[MExc].Weight >= 0 && Network[MExc-1].Weight >=0; ++MExc);
	if (MExc > 0) {
		NExc = Network[MExc-1].NStart;
	}
	else {
		NExc = 0;
	}

	// Initializing U, V
	getInputfromStruct<float>(InputmxArray, "U", this->U, getInputOps(2, "is_required", "required_size", N));
	getInputfromStruct<float>(InputmxArray, "V", this->V, getInputOps(2, "is_required", "required_size", N));

	// Initializing SpikeList
	getInputfromStruct<int32_t>(InputmxArray, "SpikeList.TimeRchd", this->SpikeList.TimeRchd,
	                            getInputOps(2, "is_required", "required_size", T));
	getInputfromStruct<int32_t>(InputmxArray, "SpikeList.TimeRchdStartInds", this->SpikeList.TimeRchdStartInds,
	                            getInputOps(2, "is_required", "required_size", T+1));

	this->StartTime = SpikeList.TimeRchd[0];
	this->StartSpikeListIndex = SpikeList.TimeRchdStartInds[0];
	this->EndSpikeListIndex = SpikeList.TimeRchdStartInds.last();

	getInputfromStruct<int32_t>(InputmxArray, "SpikeList.SpikeSynInds", this->SpikeList.SpikeSynInds,
	                            getInputOps(2, "is_required", "required_size", this->EndSpikeListIndex - this->StartSpikeListIndex));

	// Initializing Intermediate Variables
	GenSpikeList.resize(N);
	SpikeValidity.resize(N);
	SpikeTipping.resize(N);
	RespSpikeVectSplit.resize(N, FlatVectTree<uint32_t>(1));
	TippingSpikeVectSplit.resize(N, FlatVectTree<uint32_t>(1));
	CurrentGenSpikeIndex.resize(N, 0);

	// Initializing Output Variables
	ResponsibleSpikes.setDepth(2);
	TippingSpikes.setDepth(2);
	GenSpikeListOut.setDepth(1);

}

mxArray *getRespSpikes::SimulationClass::getOutput(void) {

	mxArrayPtr OutputStruct = assignmxStruct(
		{
			"ResponsibleSpikes",
			"TippingSpikes",
		    "GenSpikeList",
		},
		{
			assignmxArray(ResponsibleSpikes),
			assignmxArray(TippingSpikes),
		    assignmxArray(GenSpikeListOut),
		}
	);
	return OutputStruct;
}

void getRespSpikes::SimulationClass::ResponsibleSynCalc(void) {

	SpikeValidityCalc();
	SpikeTippingCalc();

	// CurrentRangeRespSpikes[i] is the vector of all spikes so far that contribute
	// to the spike (generated by neuron i) being iterated through. This vector
	// is used to make the process of push_back into the FLatVectArrays more efficient
	// as there is no efficient techniqu to push scalars into FlatCellArray
	MexVector<MexVector<uint32_t, CAllocator>> CurrentRangeRespSpikes(N);
	MexVector<MexVector<uint32_t, CAllocator>> CurrentRangeTippingSpikes(N);

	for(uint32_t t=StartTime; t<T+StartTime; ++t) {

		// Filtering the spikes that have arrived in the current time instant
		// according to which neuron's spike they have contributed too, or whether
		// they have contributed to any spike at all
		auto CurrentSpikeListBeg = SpikeList.TimeRchdStartInds[t-StartTime] - StartSpikeListIndex;
		auto CurrentSpikeListEnd = SpikeList.TimeRchdStartInds[t+1-StartTime] - StartSpikeListIndex;
		for(uint32_t j=CurrentSpikeListBeg; j<CurrentSpikeListEnd; ++j) {
			auto currSynapseInd = SpikeList.SpikeSynInds[j];
			auto &currSynapse = Network[currSynapseInd];
			auto currEndNeuronInd = currSynapse.NEnd-1;
			auto currentRespRange = (CurrentGenSpikeIndex[currEndNeuronInd] < SpikeValidity[currEndNeuronInd].size()) ?
				                         SpikeValidity[currEndNeuronInd][CurrentGenSpikeIndex[currEndNeuronInd]] :
				                         DiscreteRange(0,0);
			auto currentTipRange = (CurrentGenSpikeIndex[currEndNeuronInd] < SpikeTipping[currEndNeuronInd].size()) ?
				                       SpikeTipping[currEndNeuronInd][CurrentGenSpikeIndex[currEndNeuronInd]] :
				                       DiscreteRange(0,0);

			if (currentRespRange.contains(t)) {
				if (currSynapse.NStart <= NExc)
					CurrentRangeRespSpikes[currEndNeuronInd].push_back(j + StartSpikeListIndex);
			}
			if (currentTipRange.contains(t)) {
				if (currSynapse.NStart <= NExc)
					CurrentRangeTippingSpikes[currEndNeuronInd].push_back(j + StartSpikeListIndex);
			}
		}

		// Setting up CurrentGenSpikeIndex for the next iteration
		// Here, we also push elements in CurrentRangeRespSpikes into RespSpikeVectSplit
		// if the coresponding range gets over
		for(uint32_t j=0; j<N; ++j) {
			auto &GenSpikeIndex_j = CurrentGenSpikeIndex[j];
			auto &SpikeValidity_j = SpikeValidity[j];
			if (GenSpikeIndex_j < SpikeValidity[j].size()) {
				if (t == SpikeValidity_j[GenSpikeIndex_j].endPoint-1) {
					GenSpikeIndex_j++;
					RespSpikeVectSplit[j].push_back(CurrentRangeRespSpikes[j]);
					CurrentRangeRespSpikes[j].clear();
					TippingSpikeVectSplit[j].push_back(CurrentRangeTippingSpikes[j]);
					CurrentRangeTippingSpikes[j].clear();
				}
			}
		}

		if (IsProgramInterrupted()) {
			// Simply Exit. All required procedures complete in a well defined manner
			ResetInterrupt();
			break;
		}
	}

	// Here we join RespSpikeVectSplit to create ResponsibleSpikes
	for(uint32_t i=0; i<N; ++i) {
		ResponsibleSpikes.push_back(RespSpikeVectSplit[i]);
		TippingSpikes.push_back(TippingSpikeVectSplit[i]);
	}
	// We Flatten GenSpikeList to get GenSpikeListOut
	GenSpikeListOut.append(GenSpikeList);
}
