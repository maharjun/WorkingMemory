function [ IncomingSyns ] = getIncomingSyns( InputStruct )
%GETINCOMINGSYNS returns a cell array containing the incoming synapses for
%each neuron
% 
%    IncomingSyns = getIncomingSyns( InputStruct )
%   
%  Returns a cell array IncomingSyns such that IncomingSyns{i} is a column
%  vector containing the indexes of all the synapses which have post-
%  synaptic neuron as i

% Calculate Flip-Sorted Synapse Array.
[~, FSSynapses] = sortrows([InputStruct.NEnd, InputStruct.NStart]);

% Calculate the associative cell array for incoming synapses
FSNEnd = InputStruct.NEnd(FSSynapses);                                     % Getting vector of NEnds (sorted by NEnds)
[NEndExist, NEndBegLocs] = ismember((1:length(InputStruct.a))', FSNEnd);   % Find Beg of sections corrsponding to NEnd=i
NEndBegLocsExist = NEndBegLocs(NEndExist);                                 % Filtering out cases where NEnd Exists

NEndBegLocs(end+1) = length(InputStruct.NEnd);                             % Putting last element as = M
NEndBegLocsExist(end+1) = length(InputStruct.NEnd);                        

NEndBegLocsExist(1:end-1) = -diff(NEndBegLocsExist);                       % Replacing the elements except last one by
NEndBegLocs(boolean([NEndExist;1])) = NEndBegLocsExist;                    % -(no of syns belonging to that particular NEnd)

NEndBegLocs = cumsum(NEndBegLocs(end:-1:1));                               % Performing a reverse cumsum to fill in the indices
NEndBegLocs = NEndBegLocs(end:-1:1);                                       % with the actual starting index of each of the NEnd
                                                                           % Sections.
IncomingSyns = cell(length(InputStruct.a), 1);                              
for i = 1:length(InputStruct.a)                                            % Filling the Incoming synapses by partitioning
	IncomingSyns{i} = FSSynapses(NEndBegLocs(i):NEndBegLocs(i+1)-1);       % FSSynapses according to the indices in NEndBegLocs
end

end

