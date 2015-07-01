function [ A, InhNeurons, Weights, Delays ] = RecurrentNetwork(InputParams)

%RECURRENTNETWORK Summary of this function goes here
%   Detailed explanation goes here
% 
% NExc, NInh, ...
% NSynExctoExc, ...
% NSynExctoInh, ...
% NSynInhtoExc, ...
% MeanExctoExc, ...
% MeanExctoInh, ...
% MeanInhtoExc, ...
% Var, ...
% DelayRange

NExc = InputParams.NExc;
NInh = InputParams.NInh;
F_EE  = InputParams.NSynExctoExc;
F_EI  = InputParams.NSynExctoInh;
F_IE  = InputParams.NSynInhtoExc;
m_EE  = InputParams.MeanExctoExc;
m_EI  = InputParams.MeanExctoInh;
m_IE  = InputParams.MeanInhtoExc;
Var   = InputParams.Var;
DelayRange = InputParams.DelayRange;

InhNeurons = [false(NExc, 1); true(NInh, 1)];
N = InputParams.NExc + InputParams.NInh;

NEnd = cell(N, 1);
NStart = cell(N, 1);
InhSynCell = cell(N,1);
Weight = cell(N,1);
Delay = cell(N,1);

for i=1:NExc
	Temp = randperm(NExc-1, F_EE);
	Temp = [Temp(Temp < i) (Temp(Temp >= i)+1)];    % Randomly choosing w/o selfloops
	NEnd{i}       = [Temp'; randperm(NInh, F_EI)' + NExc];
	NStart{i}     = i*ones(size(NEnd{i}));
	InhSynCell{i} = false(size(NEnd{i}));
	Delay{i}      = randsample(DelayRange, F_EE + F_EI, true);
	Weight{i}     = [randn(F_EE, 1)*Var + m_EE; randn(F_EI, 1)*Var + m_EI];
	Weight{i}(Weight{i} < 0) = 0;
	if mod(i,1000) == 0
		display(i);
	end
end
for i=NExc+1:N
	NEnd{i}       = randperm(NExc, F_IE)';
	NStart{i}     = i*ones(size(NEnd{i}));
	InhSynCell{i} = true(size(NEnd{i}));
	Delay{i}      = randsample(DelayRange, F_IE, true);
	Weight{i}     = randn(F_IE, 1)*Var + m_IE;
	Weight{i}(Weight{i} > 0) = 0;
	if mod(i,1000) == 0
		display(i);
	end
end

NStartVect = cell2mat(NStart);
NEndVect = cell2mat(NEnd);
Weights = cell2mat(Weight);
Delays = cell2mat(Delay);
% clear NStart NEnd InhSynCell Weight Delay;
A = sparse(NEndVect, NStartVect, true(size(NStartVect)), N, N);

end

