function [ fig ] = PlotSpikeList( varargin )
%PLOTSPIKELIST Plots the Spike List for the given time range
%   
% There are two ways of calling this function
% 
%   1. [ fig ] = PlotSpikeList(BegTime, EndTime, InputStruct, TimeArray, SpikeList)
%   2. [ fig ] = PlotSpikeList(InputStruct, Time, SpikeList)
% 
% Return Variables:
% 
%   1. fig - Handle to generated figure
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
%   This function plots the spikes that were generated in the interval
%   [BegTime, EndTime). If not specified, the entire range of time
%   specified by TimeArray is taken into consideration.

if nargin == 5 && ~ischar(varargin{4})
	InputStruct = varargin{3};
	TimeArray   = varargin{4};
	SpikeList   = varargin{5};
	
	BegTime = varargin{1};
	EndTime = varargin{2};
elseif nargin == 3
	InputStruct = varargin{1};
	TimeArray   = varargin{2};
	SpikeList   = varargin{3};
	
	BegTime = 0;
	EndTime = (TimeArray(end) + 1)/(1000*double(InputStruct.onemsbyTstep));
end

if nargin > 5
	ExtraArgument = 6;
elseif ischar(varargin{4})
	ExtraArgument = 4;
else
	ExtraArgument = 0;
end

if ExtraArgument ~= 0 && strcmpi(varargin{ExtraArgument}, 'MarkerSize')
	MarkerSize = varargin{ExtraArgument + 1};
else
	MarkerSize = 1;
end

[GenerationTimeVect, SpikeSynIndVect] = ParseSpikeList(BegTime, EndTime, InputStruct, TimeArray, SpikeList);

SpikePreSynNeuronVect = InputStruct.NStart(SpikeSynIndVect);

fig = figure;
plot(double(GenerationTimeVect)/(1000*double(InputStruct.onemsbyTstep)), double(SpikePreSynNeuronVect), '.', 'MarkerSize', MarkerSize);
end

