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
SingleState.Iin = StateStruct.Iin(:, timeIndex);
SingleState.WeightDeriv = StateStruct.WeightDeriv(:, timeIndex);
SingleState.Weight = StateStruct.Weight(:, timeIndex);

SingleState.Iext.Iext = StateStruct.Iext.Iext(:, timeIndex);
SingleState.Iext.IExtGenState = StateStruct.Iext.IExtGenState(:, timeIndex);
SingleState.Iext.IExtNeuron = StateStruct.Iext.IExtNeuron(timeIndex);

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

