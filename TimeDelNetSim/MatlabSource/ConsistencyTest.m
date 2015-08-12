addpath('..\..\x64\Release_Exe\');


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
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.alpha = single(0.3);
InputStruct.I0 = single(1.3);
InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(80000);
InputStruct.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputStruct.StorageStepSize       = int32(1000);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);

InputStruct.OutputFile = 'SimResults1000DebugSparseLong.mat';
save('../Data/InputData.mat', 'InputStruct');

% [OutputVars, StateVars, FinalState, InitState] = TimeDelNetSim(InputStruct);
% Run the program after this
cd ..
! "..\x64\Release_Exe\TimeDelNetSim.exe"
cd MatlabSource
%% Get Detailed vector from Initial State 
% This is to check correctness of initial state return with default inputs

load('../Data/SimResults1000DebugSparseLong.mat', 'InitState');

% Setting up input settings
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
InputStruct = ConvertStatetoInitialCond(InitState);
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.alpha = single(0.3);
InputStruct.I0 = single(1.3);
InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(2000);
InputStruct.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);

InputStruct.OutputFile = 'SimResults1000DebugDetailedfromInit.mat';
save('../Data/InputData.mat', 'InputStruct');
% Run the program
cd ..
! "..\x64\Release_Exe\TimeDelNetSim.exe"
cd MatlabSource
%% Loading Relevent Data

% Loading and renaming variables for detailed simulation
load('../Data/SimResults1000DebugDetailedfromInit.mat');
clear OutputVarsDetailed StateVarsDetailed InitStateDetailed FinalStateDetailed;
OutputVarsDetailed = OutputVars;
StateVarsDetailed = StateVars;
InitStateDetailed = InitState;
FinalStateDetailed = FinalState;
clear OutputVars StateVars InitState FinalState;

% Loading and renaming variables for sparse simulation
load('../Data/SimResults1000DebugSparseLong.mat');
clear OutputVarsSparse StateVarsSparse InitStateSparse FinalStateSparse;
OutputVarsSparse = OutputVars;
StateVarsSparse = StateVars;
InitStateSparse = InitState;
FinalStateSparse = FinalState;
clear OutputVars StateVars InitState FinalState;

%% Performing Relevant Tests
max(abs(StateVarsSparse.V(:,1) - StateVarsDetailed.V(:, 4000)))

%% Getting Detailed using Final State Returned
% This is to test accurate return of final state

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
InputStruct = ConvertStatetoInitialCond(FinalStateDetailed);
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.alpha = single(0.3);
InputStruct.I0 = single(1.3);
InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(2000);
InputStruct.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);

InputStruct.OutputFile = 'SimResults1000DebugDetailedfromFinal.mat';
save('../Data/InputData.mat', 'InputStruct');
cd ..
! "..\x64\Release_Exe\TimeDelNetSim.exe"
cd MatlabSource
%% Loading Relevant Data
load('../Data/SimResults1000DebugDetailedfromFinal.mat');
clear OutputVarsDetailed StateVarsDetailed InitStateDetailed FinalStateDetailed;
OutputVarsDetailed = OutputVars;
StateVarsDetailed = StateVars;
InitStateDetailed = InitState;
FinalStateDetailed = FinalState;
clear OutputVars StateVars InitState FinalState;

%% Performing Relevant Tests
max(abs(StateVarsDetailed.V(:,4000) - StateVarsSparse.V(:,3)))

%% Getting Detailed using Intermediate Sparse State Returned
% This tests the correctness of the input of initial conditions and
% correctness of state output and state conversion to initial conditions


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
InputStruct = ConvertStatetoInitialCond(StateVarsSparse, 4*4000);
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.alpha = single(0.3);
InputStruct.I0 = single(1.3);
InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(2000);
InputStruct.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);

InputStruct.OutputFile = 'SimResults1000DebugDetailedfromInter.mat';
save('../Data/InputData.mat', 'InputStruct');
cd ..
! "..\x64\Release_Exe\TimeDelNetSim.exe"
cd MatlabSource
%% Loading Relevant Data
load('../Data/SimResults1000DebugDetailedfromInter.mat');
clear OutputVarsDetailed StateVarsDetailed InitStateDetailed FinalStateDetailed;
OutputVarsDetailed = OutputVars;
StateVarsDetailed = StateVars;
InitStateDetailed = InitState;
FinalStateDetailed = FinalState;
clear OutputVars StateVars InitState FinalState;

%% Performing Relevant Tests
max(abs(StateVarsDetailed.V(:,8000) - StateVarsSparse.V(:,6)))

%% Testing SpikeList generation

%% Get Detailed vector from Initial State 
% This is to check correctness of initial state return with default inputs

load('../Data/SimResults1000DebugSparseLong.mat', 'InitState');

% Setting up input settings
OutputOptions = { ...
	'SpikeList', ...
	'Initial'
	};

% Clearing InputStruct
clear InputStruct;

% Getting Midway state
InputStruct = ConvertStatetoInitialCond(InitState);
InputStruct.a = single(a);
InputStruct.b = single(b);
InputStruct.c = single(c);
InputStruct.d = single(d);

InputStruct.NStart = int32(NStartVect);
InputStruct.NEnd   = int32(NEndVect);
InputStruct.Weight = single(Weights);
InputStruct.Delay  = single(Delays);

InputStruct.alpha = single(0.3);
InputStruct.I0 = single(1.3);
InputStruct.onemsbyTstep          = int32(4);
InputStruct.NoOfms                = int32(80000);
InputStruct.DelayRange            = int32(RecurrentNetParams.DelayRange);
InputStruct.StorageStepSize       = int32(0);
InputStruct.OutputControl         = strjoin(OutputOptions);
InputStruct.StatusDisplayInterval = int32(8000);

InputStruct.OutputFile = 'SimResults1000DebugSpikeListfromInit.mat';
save('../Data/InputData.mat', 'InputStruct');
cd ..
! "..\x64\Release_Exe\TimeDelNetSim.exe"
cd MatlabSource
%% Loading Relevant Data

load('../Data/SimResults1000DebugSpikeListfromInit.mat');
clear OutputVarsDetailed StateVarsDetailed InitStateDetailed FinalStateDetailed;
OutputVarsDetailed = OutputVars;
StateVarsDetailed = StateVars;
InitStateDetailed = InitState;
FinalStateDetailed = FinalState;
clear OutputVars StateVars InitState FinalState;