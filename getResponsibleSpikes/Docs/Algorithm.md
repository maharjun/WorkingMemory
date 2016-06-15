Algorithm Description
=====================

##  Important Variables  ##

###   SpikeValidity   ###

Consider the DiscreteRange vectvect `SpikeValidity`

> `SpikeValidity[i][j]` is the Range of points corresponding to the jth spike
> of the ith Neuron in which an incoming spike is considered to have contribed
> to the j'th spike

###   SpikeList   ###

This is the same data structure which is used to return SpikeLists as a part 
of TimeDelNetSim. It is detailed here

    struct SpikeListStruct{
        MexVector<int32_t> SpikeSynInds;
        MexVector<int32_t> TimeRchdStartInds;
        MexVector<int32_t> TimeRchd;
    } SpikeList;

The Data in SpikeList is to be read according to the following rule:

> `SpikeSynInds[TimeRchdStartInds[j]..TimeRchdStartInds[j+1]]` (not including 
> `TimeRchdStartInds[j+1]`) represent the values of the synapses whose spike
> arrives at the time instant given by `TimeRchd[j]`.

###   CurrentGenSpikeIndex   ###

`CurrentGenSpikeIndex` is a state variable that changes with each iteration. 
For a particular iteration it is guarenteed to hold the following value at the 
behinning of ttat iteration.

> `CurrentGenSpikeIndex[i]` gives the index of the interval in SpikeValidity[i]
> whixh is either the interval that contains the current time instant or is the
> first interval that is ahead of the current time interval

###   RespSynVectSplit   ###

This is a MexVector containing "integer FlatVectTree of depth 1" i.e each 
element represents a cell array such that.

> `RespSynVectSplit[i][j]` is a vector containing the indices of all the spikes
> which were responsible for the `j`'th spike of Neuron `i`. The spike index is
> the index of the spike in SpikeList (Original SpikeList, not Truncated). This
> index is calculated using the values stored in SpikeList.TimeRchdStartIndex

The 'Split' in the name stands for the fact that the cell vector for each 
neuron is still separate and it hasn't been combined

###   RespSynFlatVect3   ###

This is the 'Joint' version of RespSynVectSplit. This is an integer FlatVectTree
of depth 3

##  Top Level Algorithm  ##

The Highest Level Flow is as:

    Calculate SpikeValidity
    for each TimeRchd:
        get CurrArrivingSpikes from SpikeList
        for each synapse in CurrArrivingSpikes:
            // take care of edge cases
            CurrNeuronRange = 
                SpikeValidity[synapse.NEnd][CurrentGenSpikeIndex[synapse.NEnd]]
            if TimeRchd is in CurrNeuronRange:
                if TimeRchd-1 is in CurrNeuronRange:
                    RespSynVectSplit[synapse.NEnd].append(synapse)
                else:
                    RespSynVectSplit[synapse.NEnd].push_back(synapse)
        for each neuron index i:
            if TimeRchd == SpikeValidity[i][CurrentGenSpikeIndex[i]].End-1:
                CurrentGenSpikeIndex[i]++;
    
    Join RespSynVectSplit to get RespSynFlatVect3.
    return RespSynFlatVect3.

##  Calculating SpikeValidity  ##

###   Input Variables   ###

    SpikeList - as defined above
    Network   - Tuple Array (NStart, NEnd, Weight, DelayinTSteps)
    V         - Matrix V(i,j) = V of jth neuron at ith time instant
    U         - Matrix U(i,j) = U of jth neuron at ith time instant

###   Intermediate Variables   ###

    GenSpikeList - MV<MV<int>> where GenSpikeList[i][j] gives time instant of 
                   the jth spike of the ith neuron

###   Output Variables   ###

    SpikeValidity - As defined Above

###   Algorithm   ###

    # Creating GenSpikeList
    for each TimeRchd in SpikeList:
        get CurrArrivingSyns;
        for each syn in CurrArrivingSyns:
            GenTime = syn.Delay - TimeRchd
            if GenTime > GenSpikeList[syn.NStart-1].last():
                GenSpikeList[syn.NStart-1].push_back(GenTime)
    
    # Creating SpikeValidity
    Initialize SpikeValidity with 0's
    for each neuron:
        for each Spike in GenSpikeList[neuron]
            search for the first time instant before the time
            instant of the spike that such that V,U are in the
            Reset Region (if none found return previous spike)
            -> CurrentRange.Begin
            CurrentRange.End = Time of Spike + 1
            SpikeValidity[neuron].push_back(CurrentRange)