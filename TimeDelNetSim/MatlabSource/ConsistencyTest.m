addpath('../../TimeDelNetSim_build/install');
% addpath('export_fig-master');
addpath ../../ExternalInputCurrent/MatlabSource/
addpath ../../MexMemoryInterfacing/MatlabSource/

%%
rng('default');
rng(25);
N = 1000;
E = 0.8;
WorkingMemNetParams.NExc = round(N*E);
WorkingMemNetParams.NInh = round(N - N*E);

WorkingMemNetParams.F_E  = 100;
WorkingMemNetParams.F_IE = 100;

WorkingMemNetParams.InitInhWeight = -5;
WorkingMemNetParams.InitExcWeight = 6;

WorkingMemNetParams.DelayRange   = 20;

[A, Ninh, Weights, Delays] = WorkingMemNet(WorkingMemNetParams);

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

%% Getting Long Sparse Vector

OutputOptions = {'FSF', 'Initial'};
% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.InitialState.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.V = single(-65*ones(N,1));
InputStruct.U = single(0.2*InputStruct.V);

InputStruct.onemsbyTstep                   = int32(1);
InputStruct.NoOfms                         = int32(80*1000);
InputStruct.DelayRange            = int32(WorkingMemNetParams.DelayRange);
InputStruct.StorageStepSize                = int32(4000);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval          = int32(2000);
InputStruct.InitialState.Iext.IExtGenState = uint32(30);

InputStruct.MaxSynWeight       = single(8);
InputStruct.ST_STDP_MaxRelativeInc = single(2.5);
InputStruct.Iext.IExtAmplitude = single(0);
InputStruct.Iext.AvgRandSpikeFreq = single(1);

if ~exist('../Data', 'dir')
    mkdir('../Data')
end
save('../Data/InputData.mat', 'InputStruct');

% [OutputVarsSparse, StateVarsSparse, FinalStateSparse, InputStateSparse] = TimeDelNetSimMEX_Lib(InputStruct);
% Run the program after this
!../../TimeDelNetSim_build/install/TimeDelNetSim ../Data/InputData.mat ../Data/SimResults1000DebugSparseLong.mat
%% Get Detailed vector from Initial State 
% This is to check correctness of initial state return with default inputs

load('../Data/SimResults1000DebugSparseLong.mat', 'InputState');

% Setting up output settings
OutputOptions = { ...
	'V', ...
	'Iin', ...
	'Itot', ...
	'Irand', ...
	'Initial', ...
	'Final'
	};

% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct = InputState;

InputStruct.NoOfms                = int32(8000);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);

InputStruct.OutputFile = 'SimResults1000DebugDetailedfromInit.mat';
save('../Data/InputData.mat', 'InputStruct');
% [OutputVarsDetailed1, StateVarsDetailed1, FinalStateDetailed1, InputStateDetailed1] = TimeDelNetSim(InputStruct);
clear functions;
% Run the program
!../../TimeDelNetSim_build/install/TimeDelNetSim ../Data/InputData.mat ../Data/SimResults1000DebugDetailedfromInit.mat
%% Loading Relevent Data

% Loading and renaming variables for detailed simulation
load('../Data/SimResults1000DebugDetailedfromInit.mat');
clear OutputVarsDetailed1 StateVarsDetailed1 InputStateDetailed1 FinalStateDetailed1;
OutputVarsDetailed1 = OutputVars;
StateVarsDetailed1 = StateVars;
InputStateDetailed1 = InputState;
FinalStateDetailed1 = FinalState;
clear OutputVars StateVars InputState FinalState;

%%
% Loading and renaming variables for sparse simulation
load('../Data/SimResults1000DebugSparseLong.mat');
clear OutputVarsSparse StateVarsSparse InputStateSparse FinalStateSparse;
OutputVarsSparse = OutputVars;
StateVarsSparse = StateVars;
InputStateSparse = InputState;
FinalStateSparse = FinalState;
clear OutputVars StateVars InputState FinalState;

%% Performing Relevant Tests
max(abs(StateVarsSparse.V(:,1) - StateVarsDetailed1.V(:, 4000)))

%% Getting Detailed using Final State Returned
% This is to test accurate return of final state

OutputOptions = { ...
	'V', ...
	'Iin', ...
	'Itot', ...
	'Final', ...
 	'Irand' ...
	'Initial', ...
	};
% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct                       = InputStateDetailed1;
InputStruct.InitialState          = FinalStateDetailed1;
InputStruct.NoOfms                = int32(8000);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);

InputStruct.OutputFile = 'SimResults1000DebugDetailedfromFinal.mat';
save('../Data/InputData.mat', 'InputStruct');

% Run Program
!../../TimeDelNetSim_build/install/TimeDelNetSim ../Data/InputData.mat ../Data/SimResults1000DebugDetailedfromFinal.mat
%% Loading Relevant Data
load('../Data/SimResults1000DebugDetailedfromFinal.mat');
clear OutputVarsDetailedFinal StateVarsDetailedFinal InputStateDetailedFinal FinalStateDetailedFinal;
OutputVarsDetailedFinal = OutputVars;
StateVarsDetailedFinal = StateVars;
InputStateDetailedFinal = InputState;
FinalStateDetailedFinal = FinalState;
clear OutputVars StateVars InputState FinalState;

%% Performing Relevant Tests
max(abs(StateVarsDetailedFinal.V(:,4000) - StateVarsSparse.V(:,3)))

%% Getting Detailed using Intermediate Sparse State Returned
% This tests the correctness of the input of initial conditions and
% correctness of state output and state conversion to initial conditions


OutputOptions = { ...
    'Itot', ...
	'V', ...
 	'U', ...
	};
% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct = InputStateSparse;
InputStruct.InitialState = getSingleState(StateVarsSparse, (12)*1000);
InputStruct.NoOfms                = int32(8000);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(2000);


% InputStruct.OutputFile = 'SimResults1000DebugDetailedfromInter.mat';
% save('../Data/InputData.mat', 'InputStruct');

[OutputVarsDetailed, StateVarsDetailed, FinalStateDetailed, InputStateDetailed] = TimeDelNetSim(InputStruct);

%% Performing Relevant Tests
max(abs(StateVarsDetailed.V(:,8000) - StateVarsSparse.V(:,5)))

%% Generating 5 Hours File.

% % The following code is commented out for safety of existing file
% InputStruct = InputStateSparse;
% 
% InputStruct.NoOfms                           = int32(5*60*60*1000);
% InputStruct.StorageStepSize                  = int32(60*1000);
% InputStruct.StatusDisplayInterval            = int32(4000);
% 
% InputStruct.MaxSynWeight                     = single(8);
% InputStateSparse.ST_STDP_EffectMaxCausal     = single(0.0);
% InputStateSparse.ST_STDP_EffectMaxAntiCausal = single(0.0);
% InputStruct.Iext.IExtAmplitude               = single(0);
% InputStruct.Iext.AvgRandSpikeFreq            = single(1);
% 
% InputStruct.OutputFile = 'SimResults1000Sparse5Hours.mat';
% 
% % Run Program
% cd ..
% ! "..\x64\Release_Exe\TimeDelNetSim.exe"
% cd MatlabSource

%% Loading 5 Hours Long Long Term STDP Simulation Result
load('../Data/SimResults1000Sparse5Hours.mat');
clear OutputVarsSparse StateVarsSparse InputStateSparse FinalStateSparse;
OutputVarsSparse = OutputVars;
StateVarsSparse = StateVars;
InputStateSparse = InputState;
FinalStateSparse = FinalState;
clear OutputVars StateVars InputState FinalState;

%% Spike Plots Generation
OutputOptions = {'PropSpikeList', 'IExt.IExtNeuron', 'Final', 'Initial'};

% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct = InputStateSparse;
InputStruct.InitialState = FinalStateSparse;
% This code is used to maintain compatibility with old data generated
% by previous code versions
if iscell(InputStruct.InitialState.SpikeQueue)
	InputStruct.InitialState.SpikeQueue = FlatCellArray.FlattenCellArray(InputStruct.InitialState.SpikeQueue);
	InputStruct.InitialState.SpikeQueue = InputStruct.InitialState.SpikeQueue.Convert2Struct();
end

% Checking initialization with different time;
InputStruct.InitialState.Time = InputStruct.InitialState.Time + 8000;

InputStruct.NoOfms                = int32(250*1000);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(2000);

InputStruct.ST_STDP_MaxRelativeInc = single(2.5);
InputStruct.Iext.IExtAmplitude = single(30);
InputStruct.Iext.AvgRandSpikeFreq = single(0.3);

% Setting IExtPattern
IExtPatternString = {
    'from 18000s to 18100s every 15s     '
    '    from 0 to 1s every 200ms          '
    '        from 100ms to 160ms generate 2'
    '        from   0ms to  60ms generate 1'
    'from 18120s onwards every 10s        '
    '    from 2s to 3s every 200ms         '
    '        from   0ms to  60ms generate 2'
    '        from 100ms to 160ms generate 1'
};
IExtPattern = getIExtPatternFromString(IExtPatternString);
IExtPattern.NeuronPatterns{end+1} = uint32([1, 60]);
IExtPattern.NeuronPatterns{end+1} = uint32([60, 1]);

InputStruct.Iext.IExtPattern = IExtPattern;

InputStruct.OutputFile = 'SimResults1000DebugSpikeListfrom5Hours.mat';
save('../Data/InputData.mat', 'InputStruct');

[OutputVarsSpikeList, StateVarsSpikeList, FinalStateSpikeList, InputStateSpikeList] = TimeDelNetSim(InputStruct);
clear functions;

%% Plotting SpikeList
BegTime = (5*60 + 0)*60 + 15;
EndTime = (5*60 + 0)*60 + 25;

figure;
[GenerationTimeVect, SpikeSynIndVect] = ParseSpikeList(BegTime, EndTime, InputStruct, OutputVarsSpikeList.PropSpikeList);
plot(GenerationTimeVect - BegTime*1000*double(InputStruct.onemsbyTstep), double(InputStruct.NStart(SpikeSynIndVect)), '.', 'MarkerSize', 1); 

%% Plotting IExt Pattern
BegTime = (5*60 + 0)*60 + 20;
EndTime = (5*60 + 0)*60 + 40;

RelInds = (StateVarsSpikeList.Time >= BegTime*InputStateSpikeList.onemsbyTstep*1000) & ...
	      (StateVarsSpikeList.Time <  EndTime*InputStateSpikeList.onemsbyTstep*1000);

figure;
plot (StateVarsSpikeList.Time(RelInds) - BegTime*InputStateSpikeList.onemsbyTstep*1000, StateVarsSpikeList.Iext.IExtNeuron(RelInds));

%% Calculating Responsible Synapses

% Performing detailed simulation of 18008s to 18016s
InputStruct = InputStateSpikeList;

InputStruct.NoOfms = int32(8000);
InputStruct.StorageStepSize = int32(0);

OutputOptions = {'U', 'V'};
InputStruct.OutputControl = strjoin(OutputOptions);

[~, StateVarsDetailed, ~] = TimeDelNetSim(InputStruct);

% Isolating relevant SpikeList
BegTime = InputStruct.InitialState.Time/1000;
EndTime = (InputStruct.InitialState.Time+InputStruct.NoOfms)/1000;
RelSpikeList = OutputVarsSpikeList.SpikeList;
BegIndex = find(RelSpikeList.TimeRchd >  BegTime*1000, 1, 'first');
EndIndex = find(RelSpikeList.TimeRchd <= EndTime*1000, 1, 'last');
RelSpikeList.TimeRchd = RelSpikeList.TimeRchd(BegIndex:EndIndex);
RelSpikeList.TimeRchdStartInds = RelSpikeList.TimeRchdStartInds(BegIndex:EndIndex+1);
RelSpikeList.SpikeSynInds = RelSpikeList.SpikeSynInds(RelSpikeList.TimeRchdStartInds(1)+1:RelSpikeList.TimeRchdStartInds(end));

clear RespSpikesInputStruct;
RespSpikesInputStruct.N = int32(length(InputStruct.a));

RespSpikesInputStruct.NStart = InputStruct.NStart;
RespSpikesInputStruct.NEnd   = InputStruct.NEnd;
RespSpikesInputStruct.Weight = InputStruct.InitialState.Weight;
RespSpikesInputStruct.Delay  = InputStruct.Delay;

RespSpikesInputStruct.onemsbyTstep = InputStruct.onemsbyTstep;
RespSpikesInputStruct.SpikeList = RelSpikeList;

RespSpikesInputStruct.U = StateVarsDetailed.U(:, 2:end);
RespSpikesInputStruct.V = StateVarsDetailed.V(:, 2:end);

InputStruct = RespSpikesInputStruct;
save('../Data/InputData.mat', 'InputStruct');
!../../TimeDelNetSim_build/install/getResponsibleSpikes ../Data/InputData.mat ../Data/OutputRespSpikes.mat

%% Loading Data
load('../Data/OutputRespSpikes.mat');

%% Running Mex Functions

RespSpikesStructMex = getResponsibleSpikes(RespSpikesInputStruct);
%% Running Tests

% Checking Consistency of result for the first spike of the first neuron
ResponsibleSpikes = FlatCellArray([], RespSpikesStruct.ResponsibleSpikes);
GenSpikeList = FlatCellArray([], RespSpikesStruct.GenSpikeList);

TestFail = MException('getResponsibleSpikes:TestFail', 'A Test has failed');


% Checking Consistency between Mex Lib and Exe (This HAS to be else
% serious shit is wrong)
if isequal(RespSpikesStructMex, RespSpikesStruct):
    fprintf('Consistency between Mex Lib and Exe Tested\n');
else
    throw(TestFail);
end

% Getting Responsible Synapses from both mex and .m and checking
RespSpikesFromMex = ResponsibleSpikes{1}{1};
RespSpikesFromFunction = getRespSpikesForSpike(StateVarsDetailed, InputStateSpikeList, OutputVarsSpikeList.SpikeList, 1, GenSpikeList{1}(1), 0);

% Testing Synapse correctness
if all(InputStateSpikeList.NEnd(OutputVarsSpikeList.SpikeList.SpikeSynInds(RespSpikesFromMex+1)+1) == 1)
    fprintf('Synapse correctness is tested for MEX Function\n');
else
    throw(TestFail);
end
if all(RespSpikesFromMex(:) == RespSpikesFromFunction(:))
    fprintf('Mex Function consistent with MATLAB Functions\n');
else
    throw(TestFail);
end
