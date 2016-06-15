#ifndef GET_RESP_SPIKES_HPP
#define GET_RESP_SPIKES_HPP

#include <MexMemoryInterfacing/Headers/MexMem.hpp>
#include <MexMemoryInterfacing/Headers/GenericMexIO.hpp>

#include <MexMemoryInterfacing/Headers/FlatVectTree/FlatVectTree.hpp>

#include <Grid2D/Headers/Range.hpp>
#include <TimeDelNetSim/Headers/Network.hpp>

namespace getRespSpikes {

class SimulationClass {
	// This class contains all the Data Members Input, Intermediate and Output
	// that are involved in a particular calculation of Responsible Spikes.
	//
	// It also contains the functions that act on this data as a part of the
	// calculation.

private:
	// Input Variables
	uint32_t onemsbyTstep;

	MexMatrix<float> U;
	MexMatrix<float> V;
	MexVector<Synapse> Network;

	struct SpikeListStruct{
		MexVector<int32_t> SpikeSynInds;
		MexVector<int32_t> TimeRchdStartInds;
		MexVector<int32_t> TimeRchd;
		SpikeListStruct() : SpikeSynInds(), TimeRchdStartInds(), TimeRchd() {}
	} SpikeList;

	// Intermediate Variables
	size_t N;
	size_t NExc;
	size_t M;
	size_t MExc;
	size_t T;
	size_t StartTime;
	size_t StartSpikeListIndex;
	size_t EndSpikeListIndex;

	MexVector<MexVector<uint32_t, CAllocator>> GenSpikeList;
	MexVector<MexVector<DiscreteRange, CAllocator>> SpikeValidity;
	MexVector<FlatVectTree<uint32_t>> RespSpikeVectSplit;
	MexVector<uint32_t> CurrentGenSpikeIndex;

	// Output Variables
	FlatVectTree<uint32_t> ResponsibleSpikes;
	FlatVectTree<uint32_t> GenSpikeListOut;

	void GenSpikeListCalc(void);
	void SpikeValidityCalc(void);

public:

	// default constructor is defined by default (All members have default constructors)

	void initialize(const mxArray *InputmxArray);
	mxArray* getOutput(void);

	void ResponsibleSynCalc(void);
};

}
#endif