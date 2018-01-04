%% Instructions:

% This program fits a simplified version of the exponential adaptive
% integrate-and-fire neuron model to electrophysiological data, as
% described in Hertaeg et al. 2012, Front. Comp. Neurosci. 

% The routine assumes that the raw data is presented in the form of a 
% structure "data" with the following compontents (see example data set
% "d_inter_071011_1"):
% * series: Recording series index (not used)
% * sess:   Recording session index
% * sweep:  Matrix with non-zero indices for each recording sweep 
%           (rows: session, columns: sweep)
% * t:      2D cell array with recording times for each sweep (in s)
% * v:      2D cell array with membrane potentials for each sweep (in V)
% 
% Each session must contain either f-I curve data (depolarizing current
% steps) or I-V curve data (hyperpolarizing current steps). Each sweep
% within a session represents data for a single current step. Such a data 
% structure is generically produced by the Patchmaster software by HEKA 
% Elektronik (we used version 2x60).

% For a data set called "X", the raw data should be contained in a file 
% called "d_<X>.mat". The function "GenerateTrainingSet.m" extracts the
% training data used for fitting and stores it in a separate file
% "Training_<X>.mat".

% The fitting routine "FitTrainingSet.m" returns:
% * initial0:   the initial values of the model parameters, estimated
%               directly from the data
% * start:      the initial values (can deviate from the initial values 
%               specified in initial0) used in the iteration steps (the number 
%               of iteration steps (total of repetitions) is specified in 
%               the parameter "numrep" see below)
% * results:    marix containing the optimized model parameters (row =
%               iteration number, column = parameters and weighted fitting 
%               error, in detail: 1=C, 2=gL, 3=EL, 4=sf, 5=Vup, 6=tcw, 
%               7=a (a=0), 8=b, 9=Vr, 10=Vth, 11=weighted fitting error)

% After fitting, the MEX file "voltage_AnalyticalModel_simpAdEx.c" can be 
% used to generate the voltage trace of the model with the specified 
% parameters in response to arbitrary input currents.


%% extract training set from raw data

clear all;
addpath('functions');

% %%%%%%%%%%% Define parameters etc. %%%%%%%%%%%%%%%%%%
folder = 'ExampleData';
FileIn = 'inter_071011_1';
FileOut = ['Training_' FileIn];
T_bounds_dep = [16 21];
T_bounds_hyper = [0.2 1];
thrs=-20;
fI_idx_all=[3 25 47];

% %%%%%%%%%%% Actual Function %%%%%%%%%%%%%%%%%%
[res,index1,index2,b] = GenerateTrainingSet(folder,FileIn,FileOut,T_bounds_dep,T_bounds_hyper,thrs,fI_idx_all);


%% fit training set
% (fitted parameters: 1=C, 2=gL, 3=EL, 4=sf, 5=Vup, 6=tcw, 7=a (a=0), 8=b, 9=Vr, 10=Vth)

clear all;
addpath('functions');

% %%%%%%%%%%% Define parameters etc. %%%%%%%%%%%%%%%%%%
file = 'ExampleData/';
filename = 'inter_071011_1';
dpt = 5;                    % number of used data points in the f-I curve (needs to be restricted if firing rates decrease at higher currents)
obs = 2;                    % index of the f-I curve to be used for fitting 
para_fix = [1,2,4,5,7];     % vector containing the model parameter numbers that are kept constant in the final optimization step 
numrep = 5;                 % number of repetitions of the optimization procedure
dev = 0.2;                  % tolerated deviation from the estimates of the initial values (20%)
plot_flag = 1;              % (1=plot fit, 0=do not show fit)
saveplot_flag = 1;          % (1=save plot ["plot_flag" has to be set to 1 as well], 0=do not save plot) 
savedata_flag = 1;          % (1=save data, 0=do not save data) 

% %%%%%%%%%%% Actual fit %%%%%%%%%%%%%%%%%%
[initial0,start,results]=FitTrainingSet(file,filename,dpt,obs,para_fix,numrep,dev,plot_flag,saveplot_flag,savedata_flag);


%% simulate model voltage trace (two examples)

% do not forget to compile the c-file with "mex":
% mex functions\voltage_AnalyticalModel_simpAdEx.c

clear all;
addpath('functions');

% define parameters and files
file='ExampleData/';                        % folder name
file_para='para_inter_071011_1_obs2_dp5';   % file with optimized parameter set
file_random='Iinj_s35_m15_out.mat';         % file with events (fluctuating current)
line=2;                                     % one of the 5 repetitions (see parameter "obs" above), e.g. 2

% load data
load([file file_para]);
load([file file_random]);

SimPar = set_parameters;
SimPar.CellPar(2,:) = results(line,1);      % C
SimPar.CellPar(3,:) = results(line,2);      % gL
SimPar.CellPar(4,:) = results(line,3);      % EL
SimPar.CellPar(5,:) = results(line,4);      % sf   
SimPar.CellPar(6,:) = results(line,5);      % Vup    
SimPar.CellPar(7,:) = results(line,6);      % tcw  
SimPar.CellPar(8,:) = results(line,7);      % a
SimPar.CellPar(9,:) = results(line,8);      % b
SimPar.CellPar(10,:) = results(line,9);     % Vr    
SimPar.CellPar(11,:) = results(line,10);    % Vth

% %%%%%%%%%%%%%%%%%%%% 1. example: constant current %%%%%%%%%%%%%%%%%%%%

Iinj=52;
Tstart=0;
Tstop=5000;

voltage_AnalyticalModel_simpAdEx(SimPar.CellPar(2:11,:),[Iinj;Tstart;Tstop],[SimPar.CellPar(4,:) 0],[0 Tstop 0.05]);
X=load('voltage_data.dat');
Tall=X(1:end,1);
Vall=X(1:end,2);

figure;
plot(Tall,Vall,'k.-');
ylim([-80 -10]);
xlim([0 Tstop]);
xlabel('t/ms'),ylabel('V/mV');

% %%%%%%%%%%%%%%%%%%%% 2. example: fluctuating current %%%%%%%%%%%%%%%%%%%%

T_bounds = [2500 Tsamp(end)-2500];
T = Tsamp(Tsamp>T_bounds(1) & Tsamp<T_bounds(2)); 
Texp = T(1:end-1);
Iapp = Iinj_samp(Tsamp>T_bounds(1) & Tsamp<T_bounds(2))*1000; Iapp = Iapp(1:end-1);
Tevents=[T_bounds(1) Texp(1:end-1)';Texp'];
events=[Iapp';Tevents];

voltage_AnalyticalModel_simpAdEx(SimPar.CellPar(2:11,:),events,[SimPar.CellPar(4,:) 0],[0 Tsamp(end) 0.05]);
X=load('voltage_data.dat');
Tall=X(1:end,1);
Vall=X(1:end,2);

figure;
plot(Tall,Vall,'k.-');
ylim([-80 -10]);
xlabel('t/ms'),ylabel('V/mV');


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
