function [ GenerationTimeVect, SpikeSynIndVect  ] = ParseSpikeList( varargin )
%PARSESPIKELIST Parses SpikeList Information and returns relevant output
%vectors
% 
% There are two ways of calling this function
% 
%   1. [ GenerationTimeVect, SpikeSynIndVect  ] = ParseSpikeList(BegTime, EndTime, InputStruct, TimeArray, SpikeList)
%   2. [ GenerationTimeVect, SpikeSynIndVect  ] = ParseSpikeList(InputStruct, Time, SpikeList)
% 
% Return Variables:
% 
%   1. GenerationTimeVect - The Vector such that GenerationTimeVect(i) is
%                           the generation time of the ith spike.
%   2. SpikeSynIndVect - The Vector such that SpikeSynIndVect(i) is the
%                        index of the synapse to which this spike belongs
% 
% Input Arguments:
% 
%   1. BegTime     - Beginning of Time Range (in seconds)
%   2. EndTime     - End of Time Range (in seconds)
%   3. InputStruct - InputStruct containing Information about the network
%                    and onemsbyTstep
%   4. TimeArray   - Array of Time Instants as returned by StateVars
%   5. SpikeList   - SpikeList Struct as returned in OutputVars.
% 
% Function Description:
% 
%   The function parses the spike list information in SpikeList, using the
%   network information in InputStruct, And returns a set of Spikes
%   Generated from BegTime to EndTime (excluding EndTime).
%   
%   The ith spike is identified by the tuple 
%   
%      (GenerationTimeVect(i), SpikeSynIndVect(i))
%   
%   The Vectors GenerationTimeVect, SpikeSynIndVect are sorted
%   lexicographically by (GenerationTimeVect(i), SpikeSynIndVect(i))
%   

if nargin == 5
	InputStruct = varargin{3};
	TimeArray   = varargin{4};
	SpikeList   = varargin{5};
	
	BegTime = varargin{1}*1000*InputStruct.onemsbyTstep;
	EndTime = varargin{2}*1000*InputStruct.onemsbyTstep;
elseif nargin == 3
	InputStruct = varargin{1};
	TimeArray   = varargin{2};
	SpikeList   = varargin{3};
	
	BegTime = TimeArray(1);
	EndTime = TimeArray(end) + 1;
end

RelTimes = TimeArray >= BegTime & TimeArray < EndTime;
BegTimeIndex = find(RelTimes, 1, 'first');
EndTimeIndex = find(RelTimes, 1, 'last') + InputStruct.onemsbyTstep*InputStruct.DelayRange;
% All Spikes generated in time corresponding to 
% 
%   TimeArray(find(RelTimes, 1, 'first'):find(RelTimes, 1, 'last'))
% 
% arrive in the time interval corresponding to
% 
%   TimeArray(find(RelTimes, 1, 'first')+1:find(RelTimes, 1, 'last') + InputStruct.onemsbyTstep*InputStruct.DelayRange);

% Calculating Total number of spikes
% TimeRchdStartInds(EndTimeIndex) = StartingIndex of Spikes landing at time
% instant TimeArray(EndTimeIndex + 1)
% Thus the above counts all spikes landing in
% TimeArray(BegTimeIndex+1:EndTimeIndex + 1)
TimeRchdStartInds = SpikeList.TimeRchdStartInds;
TotalLength = double(TimeRchdStartInds(EndTimeIndex) - TimeRchdStartInds(BegTimeIndex));

% Calculating the vector of time instants corresponding to arrival times
ArrivalTimeVect = zeros(TotalLength, 1);
InsertIndex = 1;
% By the nature of output, the entry at time t corresponds to 
% the spikes that arrive at t+1
Time = TimeArray(BegTimeIndex) + 1;
for i = BegTimeIndex:EndTimeIndex-1
	NumofElemsCurrTime = double(TimeRchdStartInds(i+1) - TimeRchdStartInds(i));
	ArrivalTimeVect(InsertIndex:InsertIndex + NumofElemsCurrTime - 1) = Time; 
	Time = Time+1;
	InsertIndex = InsertIndex + NumofElemsCurrTime;
end

% Getting Spike Information
SpikeSynIndVect       = SpikeList.SpikeSynInds(TimeRchdStartInds(BegTimeIndex)+1:TimeRchdStartInds(EndTimeIndex)) + 1;
SpikeDelayVect        = round(double(InputStruct.onemsbyTstep)*InputStruct.Delay(SpikeSynIndVect));

% Adjusting TimeVect for Delays
GenerationTimeVect    = ArrivalTimeVect - double(SpikeDelayVect);

% Filtering according to GenerationTimeVect
SpikeSynIndVect       = SpikeSynIndVect   (GenerationTimeVect >= BegTime & GenerationTimeVect <= EndTime);
GenerationTimeVect    = GenerationTimeVect(GenerationTimeVect >= BegTime & GenerationTimeVect <= EndTime);

% Sorting Spikes by (GenerationTime, SpikeSynInd)
[~, sortinds] = sortrows([GenerationTimeVect, SpikeSynIndVect]);
GenerationTimeVect    = GenerationTimeVect(sortinds);
SpikeSynIndVect       = SpikeSynIndVect(sortinds);

end

