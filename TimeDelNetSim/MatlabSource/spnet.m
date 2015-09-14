rand('seed',20);
% spnet.m: Spiking network with axonal conduction delays and STDP
% Created by Eugene M.Izhikevich.                February 3, 2004
% Modified to allow arbitrary delay distributions.  April 16,2008
M=100;                 % number of synapses per neuron
D=20;                  % maximal conduction delay 
% excitatory neurons   % inhibitory neurons      % total number 
Ne=800;                Ni=200;                   N=Ne+Ni;
a=[0.02*ones(Ne,1,'single');    0.1*ones(Ni,1,'single')];
d=[   8*ones(Ne,1,'single');    2*ones(Ni,1,'single')];
sm=10;                 % maximal synaptic strength

% post=ceil([N*rand(Ne,M);Ne*rand(Ni,M)]); 
% Take special care not to have multiple connections between neurons
delays = cell(N,D);
for i=1:Ne
    p=randperm(N-1, M);
% 	p = p(1:M);
% 	p(p>=i) = p(p>=i) + 1;
	p = [p(p<i), p(p>=i)+1];
    post(i,:)=p(1:M);
    for j=1:M
        delays{i, ceil(D*rand)}(end+1) = j;  % Assign random exc delays
    end;
end;
for i=Ne+1:N
    p=randperm(Ne);
    post(i,:)=p(1:M);
    delays{i,1}=1:M;                    % all inh delays are 1 ms.
end;

s=[6*ones(Ne,M,'single');-5*ones(Ni,M,'single')];         % synaptic weights
sd=zeros(N,M,'single');                          % their derivatives

% Make links at postsynaptic targets to the presynaptic weights
pre = cell(N,1);
aux = cell(N,1);
for i=1:Ne
    for j=1:D
        for k=1:length(delays{i,j})
            pre{post(i, delays{i, j}(k))}(end+1) = N*(delays{i, j}(k)-1)+i;
            aux{post(i, delays{i, j}(k))}(end+1) = N*(D-1-j)+i; % takes into account delay
        end;
    end;
end;
  

STDP = zeros(N,1001+D, 'single');
v = -65*ones(N,1,'single');                      % initial values
u = 0.2.*v;                             % initial values
firings=[-D 0];                         % spike timings

%%
debug;

%% Convert Network to format used in the other files.
% a = a;
b = 0.2*ones(N,1,'single');
c = -65*ones(N,1,'single');
% d = d;

NStartVect = reshape(repmat(1:N, M, 1), [], 1);
NEndVect = reshape(post', [], 1);

Delays = zeros(N, M);
for i=1:N
	for j=1:D
		for k=1:length(delays{i,j})
			Delays(i, delays{i,j}) = j;
		end
	end
end
Delays = reshape(Delays', [], 1);
Weights = reshape(s', [], 1);

