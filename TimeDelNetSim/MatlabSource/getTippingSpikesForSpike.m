function [ RespSpikes ] = getTippingSpikesForSpike(StateVars, InputStruct, SpikeList, Neuron, SpikeTime, PrevSpikeTime)
%GETTIPPINGSPIKESFORSPIKE Gets Vector of Input Spikes responsible for
%specified Spike.
%   Detailed explanation goes here

SpikeTimeTLIndex = find(StateVars.Time == SpikeTime, 1);
PrevSpikeTimeTLIndex = find(StateVars.Time == PrevSpikeTime, 1);
if isempty(PrevSpikeTimeTLIndex)
    PrevSpikeTimeTLIndex = 0;
end

% Calculating the First Reset Interval Before the current Time instant
isUnTipped = @(u, v) (25 - 4*0.04*(140-u)) > 0 && v < (-5 + sqrt(25 - 4*0.04*(140-u)))*12.5;
lv = SpikeTimeTLIndex;
lv = lv-1;
while ~isUnTipped(StateVars.U(Neuron, lv), StateVars.V(Neuron, lv)) && lv > PrevSpikeTimeTLIndex+1
    lv = lv - 1;
end
AfterPrevSpikeTime = StateVars.Time(lv);

% Calculate NExc
MExc = find(InputStruct.InitialState.Weight < 0, 1, 'first') - 1;
if isempty(MExc)
    MExc = length(InputStruct.InitialState.Weight);
end
if MExc > 0
    NExc = InputStruct.NStart(MExc-1);
else
    NExc = 0;
end

% Calculate Indices of Relevant Time Instants in SpikeList
SpikeTimeSLIndex = find(SpikeList.TimeRchd == SpikeTime, 1);
AfterPrevSpikeTimeSLIndex = find(SpikeList.TimeRchd == AfterPrevSpikeTime, 1);

RelevantSpikeInds = SpikeList.TimeRchdStartInds(AfterPrevSpikeTimeSLIndex)+1:SpikeList.TimeRchdStartInds(SpikeTimeSLIndex+1);

% Filter Indices to contain only spikes arriving at Neuron
RelevantSpikeInds = RelevantSpikeInds(InputStruct.NEnd(SpikeList.SpikeSynInds(RelevantSpikeInds)+1) == Neuron);

% Filter Indices to contain only spikes from Exc Neurons
RelevantSpikeInds = RelevantSpikeInds(InputStruct.NStart(SpikeList.SpikeSynInds(RelevantSpikeInds)+1) <= NExc);
RespSpikes = RelevantSpikeInds-1;

end