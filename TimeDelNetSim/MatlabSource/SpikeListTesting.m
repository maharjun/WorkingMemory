%% Convert Spike to spatio(neuro)-temporal Data.

BegTime = double((58)*1000*InputStruct.onemsbyTstep);
EndTime = double((60)*1000*InputStruct.onemsbyTstep);

RelTimes = StateVarsSpikeList.Time >= BegTime & StateVarsSpikeList.Time < EndTime;
BegTimeIndex = find(RelTimes, 1, 'first');
EndTimeIndex = find(RelTimes, 1, 'last') + InputStruct.onemsbyTstep*InputStruct.DelayRange - 1;

% Calculating Total number of spikes
SpikeSynInds = OutputVarsSpikeList.SpikeList.SpikeSynInds;
TimeRchdStartInds = OutputVarsSpikeList.SpikeList.TimeRchdStartInds;
TotalLength = double(TimeRchdStartInds(EndTimeIndex + 1) - TimeRchdStartInds(BegTimeIndex));

% Calculating the vector of time instants corresponding to arrival times
% (minus 1)

TimeVect = zeros(TotalLength, 1);

InsertIndex = 1;
Time = StateVarsSpikeList.Time(BegTimeIndex);
for i = BegTimeIndex:EndTimeIndex
	NumofElemsCurrTime = double(TimeRchdStartInds(i+1) - TimeRchdStartInds(i));
	TimeVect(InsertIndex:InsertIndex + NumofElemsCurrTime - 1) = Time;
	Time = Time+1;
	InsertIndex = InsertIndex + NumofElemsCurrTime;
end

% Straightening out the Cell Array.
SpikeListVect = SpikeSynInds(TimeRchdStartInds(BegTimeIndex)+1:TimeRchdStartInds(EndTimeIndex+1)) + 1; % +1 for the C++ to matlab 
                                         % indexing convention conversion

% Calculating Synapse parameter vectors
SpikePreSynNeuronVect = InputStruct.NStart(SpikeListVect);
SpikeDelayVect        = round(double(InputStruct.onemsbyTstep)*InputStruct.Delay(SpikeListVect));

% Adjusting TimeVect for Delays
TimeVect = TimeVect - SpikeDelayVect + 1;

% Removing all entries that do not come into the relevant time frame
SpikeListVect         = SpikeListVect(TimeVect >= BegTime & TimeVect <= EndTime);
SpikePreSynNeuronVect = SpikePreSynNeuronVect(TimeVect >= BegTime & TimeVect <= EndTime);
SpikeDelayVect        = SpikeDelayVect(TimeVect >= BegTime & TimeVect <= EndTime);
TimeVect              = TimeVect(TimeVect >= BegTime & TimeVect <= EndTime);

% Creating Sparse Matrix
TimeRange = EndTimeIndex - BegTimeIndex + 1;
SpikeMat = sparse(double(SpikePreSynNeuronVect), double(TimeVect) - BegTime + 1, ones(length(TimeVect), 1), N, double(TimeRange));
figure;

plot(double(TimeVect - min(TimeVect)), double(SpikePreSynNeuronVect), '.', 'MarkerSize', 1);

% spy(SpikeMat, 5);

%% Random Plotting
RelNeuron = 810;
PlotTime = TimeVect(SpikePreSynNeuronVect == RelNeuron) - min(TimeVect);
PlotSpikes = ones(length(PlotTime), 1);

plot(PlotTime, PlotSpikes, '.', 'MarkerSize', 5);

