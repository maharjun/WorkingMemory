function [ EffectiveWeights ] = getEffectiveWeights(Weights, ST_STDP_RelativeInc)
%GETEFFECTIVEWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

EffectiveWeights = ST_STDP_RelativeInc;
EffectiveWeights(EffectiveWeights < -1) = -1;
EffectiveWeights = Weights.*(1+EffectiveWeights);

end

