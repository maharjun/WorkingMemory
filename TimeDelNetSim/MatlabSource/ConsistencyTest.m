rmpath('../../x64/Debug_Lib');
addpath('../../x64/Release_Lib');
% addpath('export_fig-master');
addpath ../Headers/IExtHeaders/MatlabSource/
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

InputStruct.OutputFile = 'SimResults1000DebugSparseLong.mat';
save('../Data/InputData.mat', 'InputStruct');

% [OutputVarsSparse, StateVarsSparse, FinalStateSparse, InputStateSparse] = TimeDelNetSimMEX_Lib(InputStruct);
% Run the program after this
! start "TimeDelNetSim Sparse Simulation" /d . "powershell" ". .\Release_Exe.ps1"
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
! start "TimeDelNetSim Sparse Simulation" /d . "powershell" ". .\Release_Exe.ps1"
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
! start "TimeDelNetSim Sparse Simulation" /d . "powershell" ". .\Release_Exe.ps1"
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
IExtPattern = getEmptyIExtPattern();

IExtPattern.NeuronPatterns{end+1} = uint32([1, 60]);
IExtPattern.NeuronPatterns{end+1} = uint32([60, 1]);

IExtPattern = AddInterval(IExtPattern, 0, 18000000, 18000000+100000, 15000, 0);      % 1
	IExtPattern = AddInterval(IExtPattern, 1, 0, 1000, 200, 0);      % 2
		IExtPattern = AddInterval(IExtPattern, 2, 100, 160 , 0, 2);   % 3
		IExtPattern = AddInterval(IExtPattern, 2,   0,  60, 0, 1);   % 4
IExtPattern = AddInterval(IExtPattern, 0, 18000000+120000, 18000000+120000, 10000, 0); % 5
	IExtPattern = AddInterval(IExtPattern, 5, 2000, 3000, 200, 0);   % 6
		IExtPattern = AddInterval(IExtPattern, 6, 0  , 60 , 0, 2);   % 7
		IExtPattern = AddInterval(IExtPattern, 6, 100, 160, 0, 1);   % 8

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

plot (StateVarsSpikeList.Time(RelInds) - BegTime*InputStateSpikeList.onemsbyTstep*1000, StateVarsSpikeList.Iext.IExtNeuron(RelInds));