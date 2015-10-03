function [ OutgoingSyns ] = getOutgoingSyns( InputStruct )
%GETOUTGOINGSYNS returns a cell array containing the incoming synapses for
%each neuron
% 
%    OutgoingSyns = GetIncomingSyns( InputStruct )
%   
%  Returns a cell array OutgoingSyns such that OutgoingSyns{i} is a column
%  vector containing the indexes of all the synapses which have pre-
%  synaptic neuron as i

% Calculate the associative cell array for outgoing synapses
NStart = InputStruct.NStart;                                                 % Getting vector of NStart (sorted by NStart)
[NStartExist, NStartBegLocs] = ismember((1:length(InputStruct.a))', NStart); % Find Beg of sections corrsponding to NStart=i
NStartBegLocsExist = NStartBegLocs(NStartExist);                             % Filtering out cases where NStart=i Exists

NStartBegLocs(end+1) = length(InputStruct.NStart);                           % Putting last element as = M
NStartBegLocsExist(end+1) = length(InputStruct.NStart);                      

NStartBegLocsExist(1:end-1) = -diff(NStartBegLocsExist);                     % Replacing the elements except last one by
NStartBegLocs(boolean([NStartExist;1])) = NStartBegLocsExist;                % -(no of syns belonging to that particular NStart)

NStartBegLocs = cumsum(NStartBegLocs(end:-1:1));                             % Performing a reverse cumsum to fill in the indices
NStartBegLocs = NStartBegLocs(end:-1:1);                                     % with the actual starting index of each of the NEnd
                                                                             % Sections.
OutgoingSyns = cell(length(InputStruct.a), 1);
for i = 1:length(InputStruct.a)                                              % Filling the Incoming synapses by partitioning
	OutgoingSyns{i} = (NStartBegLocs(i):NStartBegLocs(i+1)-1)';              % 1:M according to the indices in NStartBegLocs
end

end

