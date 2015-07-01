function [ InputStruct ] = ConvertStatetoInitialCond( StateStruct, timeInstant )
% CONVERTSTATETOINITIALCOND Converts structs
%   basically a name conversion Function as such. 

timeDimensionLength = length(StateStruct.Time);
timeIndex = 1;
if timeDimensionLength == 1
	timeIndex = 1;
elseif nargin() == 2 
	timeIndex = find(StateStruct.Time == timeInstant, 1);
	if isempty(timeIndex)
		Ex = MException('NeuralSim:ConvertStatetoInitialCond:InvalidTimeInstant', ...
					'Need to specify a time instant belonging to from StateStruct.Time');
		throw(Ex);
	end

end

InputStruct.V = StateStruct.V(:, timeIndex);
InputStruct.U = StateStruct.U(:, timeIndex);
InputStruct.Iin1 = StateStruct.Iin1(:, timeIndex);
InputStruct.Iin2 = StateStruct.Iin2(:, timeIndex);
InputStruct.Irand = StateStruct.Irand(:, timeIndex);
InputStruct.GenState = StateStruct.GenState(:, timeIndex);
InputStruct.Time = StateStruct.Time(timeIndex);

InputStruct.CurrentQIndex = StateStruct.CurrentQIndex(timeIndex);
if timeDimensionLength == 1
	InputStruct.SpikeQueue = StateStruct.SpikeQueue;
else
	InputStruct.SpikeQueue = StateStruct.SpikeQueue{timeIndex};
end

InputStruct.LSTNeuron = StateStruct.LSTNeuron(:, timeIndex);
InputStruct.LSTSyn = StateStruct.LSTSyn(:, timeIndex);

end

