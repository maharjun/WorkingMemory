%%
rmpath('F:\Users\Arjun\Desktop\Acads\SRE\TimeDelNetSimMEX\x64\Debug_Lib');
addpath('F:\Users\Arjun\Desktop\Acads\SRE\TimeDelNetSimMEX\x64\Release_Lib');
addpath('export_fig-master');

%%
rng(25);
N = 1000;
E = 0.8;
RecurrentNetParams.NExc = round(N*E);
RecurrentNetParams.NInh = round(N - N*E);

RecurrentNetParams.NSynExctoExc = ceil(100*N/2000);
RecurrentNetParams.NSynExctoInh = ceil(100*N/2000);
RecurrentNetParams.NSynInhtoExc = ceil(1200*N/2000);

RecurrentNetParams.MeanExctoExc = 0.5*2000/N;
RecurrentNetParams.MeanExctoInh = 0.15*2000/N;
RecurrentNetParams.MeanInhtoExc = -0.7*2000/N;

RecurrentNetParams.Var          = 0.2;
RecurrentNetParams.DelayRange   = 20;

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
% OutputOptions = { ...
% 	'V', ...
% 	'Iin', ...
% 	'Itot', ...
% 	'Irand', ...
% 	'Initial'
% 	};
OutputOptions = {'FSF', 'Initial'};
% OutputOptions = {'SpikeList', 'Final', 'Initial'};
% Clearing InputStruct
clear InputStruct;

% Getting Midway state
% InputStruct = ConvertStatetoInitialCond(FinalState);
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(80000);
InputStruct.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputStruct.StorageStepSize       = int32(1000);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);

% tic;
% try 
% 	[OutputVar2, StateVars2, FinalState2, InitState2] = TimeDelNetSimMEX_Lib(InputStruct);
% catch e
% 	clear functions;
% 	throw(e);
% end
% toc;
InputStruct.OutputFile = 'SimResults1000DebugSparseLong.mat';
save('../Data/InputData.mat', 'InputStruct');
clear functions;
% %% Consistency Check
% ConsistencyCheck;

%% Data Load

load('TimeDelNetSimMEX\TimeDelNetSimMEX_Exe\Data\SimResults1000DebugDetailLong.mat');
% % 
% clear OutputVars2;
% OutputVars2 = OutputVars;
% clear OutputVars;
% 
% clear StateVars2;
% StateVars2 = StateVars;
% clear StateVars;
% 
% clear InitState2;
% InitState2 = InitState;
% clear InitState;

%% Grid plot
relinds = 1:InputStruct.onemsbyTstep*1000;
IClipped = OutputVarsDetailed.Itot(:,relinds);
IThresh = 30;
IClipped(IClipped > IThresh) = IThresh;
IClipped(IClipped < 0) = 0;

VClipped = StateVarsDetailed.V;
VClipped = VClipped + 60;
VClipped(VClipped > 90) = 90;
VClipped(VClipped < 0) = 0;

% StateVars2.V(:,relinds) ~= 30
% plotMat = double(VClipped/90);
plotMat = double(IClipped/IThresh);		% Black for spike
% plotMat = double(StateVars2.V(:,relinds) ~= 30);
plotMat(1:1) = 0;
[Rows,Cols] = size(plotMat);                           %# Get the matrix size

ExpandFactor = 1;
plotMatExp = zeros(size(plotMat,1)*ExpandFactor, size(plotMat,2));
for i=1:size(plotMat,1)
	plotMatExp((i-1)*ExpandFactor+1:i*ExpandFactor, :) = repmat(plotMat(i, :), ExpandFactor, 1);
end
DefFigColor = get(0, 'DefaultFigureColor');
set(0, 'DefaultFigureColor', 'w');
FigHandle = DisplayImage(plotMatExp, 'axes', true, 'Tag', 'export_fig_native');
colormap(flipud(colormap));
XTicks = get(gca, 'XTick');
% YTicks = get(gca, 'YTick');
% set(gca, 'XTickLabels', strsplit(num2str(XTicks/single(InputStruct.onemsbyTstep))), ...
% 		 'YTickLabels', strsplit(num2str(YTicks/ExpandFactor)));
C = colorbar();
CTicks = get(C, 'YTick');
set(C, 'YTickLabels', strsplit(num2str(CTicks*IThresh)));
% axis equal                                   %# Make axes grid sizes equal
% set(gca, ...
%         'GridLineStyle','-','XGrid','off','YGrid','off');
export_fig('PeriodicSpiking.png', FigHandle, '-native', '-q101', '-a1');

%% Random plots
relInds = 1:InputStruct.onemsbyTstep*InputStruct.NoOfms;
% PlotVect = StateVars.V(34, relInds);
PlotVect = OutputVarsDetailed.Iin(34, relInds);
figure; plot(StateVarsDetailed.Time(relInds), PlotVect );
mean((OutputVarsDetailed.Itot(400, relInds) - StateVarsDetailed.Iin(400, relInds)).^2)
mean((PlotVect-mean(PlotVect)).^2)