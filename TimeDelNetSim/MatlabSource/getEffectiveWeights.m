function [ EffectiveWeights ] = getEffectiveWeights(State, InputStruct)
%GETEFFECTIVEWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

EffectiveWeights = State.ST_STDP_RelativeInc;

% Calculating actual effective weights
LastUpdatedTime = max(State.LSTSyn, State.LSTNeuron(InputStruct.NEnd));
CurrTime = State.Time + 1; % we're looking at weight at end of iteration
EffectiveWeights(LastUpdatedTime > -1) = EffectiveWeights(LastUpdatedTime > -1).*(InputStruct.ST_STDP_DecayWithTime.^single(CurrTime - LastUpdatedTime(LastUpdatedTime > -1)));

EffectiveWeights(EffectiveWeights < -1) = -1;
EffectiveWeights = State.Weight.*(1+EffectiveWeights);

end

