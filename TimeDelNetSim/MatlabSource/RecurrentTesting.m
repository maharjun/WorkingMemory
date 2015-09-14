%%
if strcmpi(Build, 'Debug')
	rmpath('../../x64/Release_Lib');
	addpath('../../x64/Debug_Lib');
elseif strcmpi(Build, 'Release')
	rmpath('../../x64/Debug_Lib');
	addpath('../../x64/Release_Lib');
end
%%
rng(25);
N = 2000;
E = 0.8;
RecurrentNetParams.NExc = round(N*E);
RecurrentNetParams.NInh = round(N - N*E);

RecurrentNetParams.NSynExctoExc = 1000;
RecurrentNetParams.NSynExctoInh = 300;
RecurrentNetParams.NSynInhtoExc = 500;

RecurrentNetParams.MeanExctoExc = 0.1;
RecurrentNetParams.MeanExctoInh = 0.09;
RecurrentNetParams.MeanInhtoExc = -1;

RecurrentNetParams.Var          = 0.001;
RecurrentNetParams.DelayRange   = 30;

[A, Ninh, Weights, Delays] = RecurrentNetwork(RecurrentNetParams);

a = 0.02*ones(N,1);
b = 0.2*ones(N,1);
c = -65*ones(N,1);
d = 8*ones(N,1);

a(Ninh) = 0.1;
b(Ninh) = 0.2;
c(Ninh) = -65;
d(Ninh) = 2;
% Delays = Delays + 10;
[NEndVect, NStartVect] = find(A);

%% Input setup
% Setting up input settings
OutputOptions = {
	'V', ...
	'Iin1', ...
	'Iin2', ...
	'Final' ...
};

% Clearing InputList
clear InputList;

% Getting Midway state
% InputList = ConvertStatetoInitialCond(StateVars1, 1);
InputList.a = single(a);
InputList.b = single(b);
InputList.c = single(c);
InputList.d = single(d);

InputList.NStart = int32(NStartVect);
InputList.NEnd   = int32(NEndVect);
InputList.Weight = single(Weights);
InputList.Delay  = single(Delays);

InputList.onemsbyTstep          = int32(4);
InputList.NoOfms                = int32(200);
InputList.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputList.StorageStepSize       = int32(0);
InputList.OutputControl         = strjoin(OutputOptions);
InputList.StatusDisplayInterval = int32(600);

tic;
try 
	[OutputVar2, StateVars2, FinalState2, InitState2] = TimeDelNetSimMEX_Lib(InputList);
catch e
	clear functions;
	throw(e);
end
toc;
clear functions;

% %% Grid plot
% relinds = 1:InputList.onemsbyTstep*InputList.NoOfms;
% IClipped = StateVars2.I;
% IClipped(IClipped > 30) = 30;
% IClipped(IClipped < 0) = 0;
% VClipped = StateVars2.V;
% VClipped = VClipped + 60;
% VClipped(VClipped > 90) = 90;
% VClipped(VClipped < 0) = 0;

% % StateVars2.V(:,relinds) ~= 30
% % plotMat = double(VClipped/90);
% plotMat = double(IClipped/30);		% Black for spike
% % plotMat = double(StateVars2.V(:,relinds) ~= 30);
% plotMat(1:1) = 0;
% [Rows,Cols] = size(plotMat);                           %# Get the matrix size

% DefFigColor = get(0, 'DefaultFigureColor');
% set(0, 'DefaultFigureColor', 'w');
% FigHandle = DisplayImage(plotMat, 'axes', true, 'Tag', 'export_fig_native', 'cbar', false);
% % axis equal                                   %# Make axes grid sizes equal
% % set(gca, ...
% %         'GridLineStyle','-','XGrid','off','YGrid','off');
% export_fig('PeriodicSpiking.png', FigHandle, '-native', '-q101', '-a1');
% %% Random plots
% relInds = 1:InputList.onemsbyTstep*InputList.NoOfms;
% figure; plot(StateVars2.Time(relInds), StateVars2.U(1000, relInds));