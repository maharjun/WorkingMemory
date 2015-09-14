%%
IExt = OutputVarsDetailed.IExtNeuron;
ItotIzhik = zeros(N, 1000, 'single');
VIzhik = zeros(N, 1000, 'single');
UIzhik = zeros(N, 1000, 'single');
WeightDerivIzhik = zeros(N*M, 1000, 'single');

fired = [];
figure;
for sec=1:1                      % simulation of 1 day
  for t=1:1000                          % simulation of 1 sec
    I=zeros(N,1, 'single');        
    I(IExt((sec-1)*1000 + t))=single(20);                 % random thalamic input 
    fired = find(v>=30);                % indices of fired neurons
    v(fired)=-65;  
    u(fired)=u(fired)+d(fired);
    STDP(fired,t+D)=single(0.1);
    for k=1:length(fired)
      sd(pre{fired(k)})=sd(pre{fired(k)})+STDP(N*t+aux{fired(k)});
    end;
    firings=[firings;t*ones(length(fired),1),fired];
    k=size(firings,1);
	
	WeightDerivIzhik(:,t) = reshape(sd', [], 1);
    while firings(k,1)>t-D
      del=delays{firings(k,2),t-firings(k,1)+1};
      ind = post(firings(k,2),del);
      I(ind)=I(ind)+s(firings(k,2), del)';
      sd(firings(k,2),del)=sd(firings(k,2),del)-single(1.2)*STDP(ind,t+D)';
      k=k-1;
    end;
	ItotIzhik(:,t) = I;
	VIzhik(:,t) = v;
 	UIzhik(:,t) = u;
	
    v=v+single(0.5)*((single(0.04)*v+single(5)).*v+single(140)-u+I);    % for numerical 
	v(v < -100) = -100;
    v=v+single(0.5)*((single(0.04)*v+single(5)).*v+single(140)-u+I);    % stability time 
	v(v < -100) = -100;
    u=u+a.*(single(0.2)*v-u);               % step is 0.5 ms
 	
% 	v_new=v+((single(0.04)*v+5).*v+single(140)-u+I);
% 	u_new=u+a.*(single(0.2)*v-u);
% 	v_new(v_new < -100) = -100;
% 	v = v_new;
% 	u = u_new;
% 	v(fired)=-65;  
%     u(fired)=u(fired)+d(fired);
% 	
% 	fired = find(v>=30);
% 	v(fired)=30;  
    
    STDP(:,t+D+1)=single(0.95)*STDP(:,t+D);     % tau = 20 ms
  end;
  plot(firings(:,1),firings(:,2),'.','MarkerSize',1);
  title(sprintf('Time = %d', sec));
  axis([0 1000 0 N]); drawnow;
  STDP(:,1:D+1)=STDP(:,1001:1001+D);
  ind = find(firings(:,1) > 1001-D);
  firings=[-D 0;firings(ind,1)-1000,firings(ind,2)];
  s(1:Ne,:)=max(0,min(sm,0.01+s(1:Ne,:)+sd(1:Ne,:)));
  sd=0.9*sd;
end;