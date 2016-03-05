function [ SingleState ] = getSingleState( StateStruct, timeInstant, NoWarnings)
% CONVERTSTATETOINITIALCOND Converts structs
%   basically a name conversion Function as such. 

if nargin() == 2
	NoWarnings = false;
end

timeIndex = find(StateStruct.Time == timeInstant, 1);
if isempty(timeIndex)
	Ex = MException('NeuralSim:ConvertStatetoInitialCond:InvalidTimeInstant', ...
				'Need to specify a time instant belonging to from StateStruct.Time');
	throw(Ex);
end

SingleState = getSingleRecord(StateStruct, timeIndex, NoWarnings);

end

