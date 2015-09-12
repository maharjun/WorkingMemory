function [ SingleState ] = getSingleState( StateStruct, timeInstant )
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

SingleState.V = StateStruct.V(:, timeIndex);
SingleState.U = StateStruct.U(:, timeIndex);
SingleState.Iin1 = StateStruct.Iin1(:, timeIndex);
SingleState.Iin2 = StateStruct.Iin2(:, timeIndex);
SingleState.Iext.Irand = StateStruct.Iext.Irand(:, timeIndex);
SingleState.Iext.GenState = StateStruct.Iext.GenState(:, timeIndex);
SingleState.Time = StateStruct.Time(timeIndex);

SingleState.CurrentQIndex = StateStruct.CurrentQIndex(timeIndex);
if timeDimensionLength == 1
	SingleState.SpikeQueue = StateStruct.SpikeQueue;
else
	SingleState.SpikeQueue = StateStruct.SpikeQueue{timeIndex};
end

SingleState.LSTNeuron = StateStruct.LSTNeuron(:, timeIndex);
SingleState.LSTSyn = StateStruct.LSTSyn(:, timeIndex);

end

