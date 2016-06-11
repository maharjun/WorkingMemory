function [ SparseSpikeMat, TimeInstantVect] = getSparseSpikeArr(InputStruct, SpikeList, BegTime, EndTime)
%GETSPARSESPIKEARR Parses SpikeList Information and returns the sparse
%matrix containing spikes
% 
% There are two ways of calling this function
% 
%   1. [ SparseSpikeMat, TimeInstantVect] = ParseSpikeList(InputStruct, SpikeList)
%   2. [ SparseSpikeMat, TimeInstantVect] = ParseSpikeList(InputStruct, SpikeList, BegTime, EndTime)
% 
% Return Variables:
% 
%   1. SparseSpikeMat  - A boolean sparse matrix such that SparseSpikeMat
%                        (n,m) = 1 or 0 depending on whether it corresponds
%                        to a spike or not
%   2. TimeInstantVect - The Vector such that TimeInstantVect(i) is the 
%                        time instant corresponding to SparseSpikeMat(:,i)
% 
% Input Arguments:
% 
%   1. BegTime     - Beginning of Time Range (in seconds)
%   2. EndTime     - End of Time Range (in seconds)
%   3. InputStruct - InputStruct containing Information about the network
%                    and onemsbyTstep
%   5. SpikeList   - SpikeList Struct as returned in OutputVars.
% 
% Function Description:
% 
%   The function parses the spike list information in SpikeList, using the
%   network information in InputStruct, And returns a sparse boolean matrix
%   corresponding to the spiking activity of the given time interval
%   
%   Each entry (n,m) in SparseSpikeMat corresponds to the Neuron 'n' and a
%   generation time of TimeInstantVect(m)
%   

if nargin == 2
	BegTime = double(SpikeList.TimeRchd(1));
	EndTime = double(SpikeList.TimeRchd(end) + 1);
elseif nargin == 4
	BegTime = round(BegTime*double(InputStruct.onemsbyTstep*1000));
	EndTime = round(EndTime*double(InputStruct.onemsbyTstep*1000));
end

TimeScaleFactor = double(InputStruct.onemsbyTstep*1000);
[GenerationTimeVect, SpikeSynIndVect] = ParseSpikeList(BegTime/TimeScaleFactor, EndTime/TimeScaleFactor, InputStruct, SpikeList);

NNeurons = double(max(max(InputStruct.NStart), max(InputStruct.NEnd)));

SparseSpikeMat = sparse(...
	double(InputStruct.NStart(SpikeSynIndVect)),  ...
	GenerationTimeVect-BegTime+1, ...
	ones(length(GenerationTimeVect),1), ...
	NNeurons, EndTime-BegTime);

TimeInstantVect = (BegTime:EndTime-1)';

end

