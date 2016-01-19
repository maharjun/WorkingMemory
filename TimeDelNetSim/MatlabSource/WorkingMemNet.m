function [ A, InhNeurons, Weights, Delays ] = WorkingMemNet( NetworkParams )
%WORKINGMEMNET Summary of this function goes here
%   Detailed explanation goes here

NExc = NetworkParams.NExc; % 800
NInh = NetworkParams.NInh; % 200

InhNeurons = [false(NExc, 1); true(NInh, 1)];
N = NExc + NInh;

DelayRange = NetworkParams.DelayRange; % 20

F_E  = NetworkParams.F_E;  % Fanout of Exc Neurons = 100
F_IE = NetworkParams.F_IE; % Fanout of Inh->Exc connection = 100

InitInhWeight = NetworkParams.InitInhWeight;
InitExcWeight = NetworkParams.InitExcWeight;

NEnd = cell(N, 1);
NStart = cell(N, 1);
InhSynCell = cell(N,1);
Weight = cell(N,1);
Delay = cell(N,1);
for i=1:NExc
	Temp = randperm(N-1, F_E);
	Temp = [Temp(Temp < i) (Temp(Temp >= i)+1)];    % Randomly choosing w/o selfloops
	NEnd{i}       = Temp';
	NStart{i}     = i*ones(size(NEnd{i}));
	InhSynCell{i} = false(size(NEnd{i}));
	Delay{i}      = randsample(DelayRange, F_E, true);
	Weight{i}     = InitExcWeight*ones(size(NEnd{i}));
	Weight{i}(Weight{i} < 0) = 0;
	if mod(i,1000) == 0
		display(i);
	end
end
for i=NExc+1:N
	NEnd{i}       = randperm(NExc, F_IE)';
	NStart{i}     = i*ones(size(NEnd{i}));
	InhSynCell{i} = true(size(NEnd{i}));
	Delay{i}      = 1*ones(F_IE, 1);
	Weight{i}     = InitInhWeight*ones(size(NEnd{i}));
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

