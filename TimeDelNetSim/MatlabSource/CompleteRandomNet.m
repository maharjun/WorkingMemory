function [A, InhSyn, NonInhSyn] = CompleteRandomNet(N, Ninh)
%COMPLETERANDOMNET Summary of this function goes here
%   Detailed explanation goes here
NEnd = cell(N, 1);
NStart = cell(N, 1);
InhSynCell = cell(N,1);
for i=1:N
	Atemp = ~logical(floor(rand(N-i, 1)*50));
	NEnd{i} = find(Atemp)+i;
	NStart{i} = i*ones(size(NEnd{i}));
	if Ninh(i)
		InhSynCell{i} = true(size(NEnd{i}));
	else
		InhSynCell{i} = false(size(NEnd{i}));
	end
	if mod(i,1000) == 0
		display(i);
	end
end
NStartVect = cell2mat(NStart);
NEndVect = cell2mat(NEnd);
InhSynVect = cell2mat(InhSynCell);
clear NStart NEnd InhSynCell;
A = sparse(NEndVect, NStartVect, true(size(NStartVect)), N, N);
InhSyn = sparse(NEndVect, NStartVect, InhSynVect, N, N);
NonInhSyn = sparse(NEndVect, NStartVect, ~InhSynVect, N, N);

end

