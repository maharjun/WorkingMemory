function [ SingleRecord ] = getSingleRecord( InputData, Index )
%GETSINGLERECORD Summary of this function goes here
%   Detailed explanation goes here

if isstruct(InputData) && isfield(InputData, 'ClassName') && strcmp(InputData.ClassName, 'FlatCellArray')
	tempFCA = FlatCellArray([], InputData);
	SingleRecord = tempFCA{Index}.Convert2Struct();
elseif isstruct(InputData)
	% Applying the same function on each of its fields
	DataFieldNames = fieldnames(InputData);
	for i = 1:length(DataFieldNames)
		CurrentFieldName = DataFieldNames{i};
		SingleRecord.(CurrentFieldName) = getSingleRecord(InputData.(CurrentFieldName), Index);
	end
elseif iscell(InputData)
	% Selecting element of cell array
	SingleRecord = InputData{Index};
elseif min(size(InputData)) > 1
	% Selecting column in matrix
	SingleRecord = InputData(:, Index);
elseif isvector(InputData)
	% Selecting element of vector
	SingleRecord = InputData(Index);
end


end

