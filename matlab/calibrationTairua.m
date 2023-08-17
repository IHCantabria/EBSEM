%% Tairua calibration

% Models Application for Tairúa data

clc; clear *;
addpath ./data/
addpath ./functions/
addpath ./functions/Calibration
addpath ./functions/Waves/

load("Wave_hindcast_corrected.mat")
load("hindcast_SS_corr.mat")
load("forecast_SS.mat"); timesf = time;
load("Shorecast.mat")
load("Tide_past.mat")

%% Setting calibration period

dt = 3;
time = (Shorecast.time(1):dt/24:Shorecast.time(end))';
AT = interp1(Tide_past.time,Tide_past.tide,time);
SS = interp1([Storm_surge.time; timesf], [Storm_surge.SS; Storm_Surge], time);
Hs = interp1(hindcast.time,hindcast.Hs,time);
Tp = interp1(hindcast.time,hindcast.Tp,time);
theta = interp1(hindcast.time,hindcast.Dir,time);

ENS.Yobs = Shorecast.average;
ENS.time = Shorecast.time;

%% Params setup

d50 = 0.3e-3;
Hberm = 1;
Yi = Shorecast.average(1);
flagP = 4;
depth = 10;
angleBathy = 54.3;


for i = 1:length(ENS.time)
    [~, ii(i)] =  min(abs(ENS.time(i)-time));
end
ENS.indexes = ii;




%% Initialize class for run the models

Setup = EBSM();
Setup = Setup.init(time, Hs, SS, AT, Tp, d50);

Setup = Setup.LinearBreak(theta, depth, angleBathy);

Setup = Setup.MillerDean(dt, Yi, Hberm, flagP);
Setup = Setup.Yates09(dt, Yi);
Setup = Setup.ShoreFor(dt, Yi);

ngen = 2000;
npop = 100;
nobj = 3;
mag = 0.25;


%% Miller and Dean 2004


disp("Starting Miller and Dean 2004...")
obj = @(X) Objective("MD", Setup, X, "RMSE", ENS);
% obj = @(X) MultiObjective("MD", Setup, X, ENS);
tic;

npar = 3;
x0 = log([78.75 6.9211e-05 8.4790e-04]);
lb = log([60 5e-6 5e-5]);
ub = log([100 1e-3 8e-3]);

% [population, objectives] = nsga2(obj, npar, ngen, npop, ub, lb, nobj);
% 
% [~, bestSol] = min(sum(objectives,2));
% 
% params_md = population(bestSol,:);
% met_md = objectives(bestSol,:);

[params_md, met_md] = sce_ua2(obj, x0,ngen, npop, npar, mag, lb, ub);

disp(strcat('Elapsed time: ', string(toc), 's'));

[Y_md]= millerDean04(Setup, exp(params_md));
%% Yates 2009
% 
npar = 4;

x0 = log([0.1143 9.6392 0.0034 0.0038]);
lb = log([0.01 1 1e-4 0.001]);
ub = log([0.3 15 0.05 0.01]);

disp("Starting Yates et al. 2009...")
obj = @(X) Objective("Y09", Setup, X, "RMSE", ENS);
% obj = @(X) MultiObjective("MD", Setup, X, ENS);
tic;


% [population, objectives] = nsga2(obj, npar, ngen, npop, ub, lb, nobj);
% 
% [~, bestSol] = min(sum(objectives,2));
% 
% params_y09 = population(bestSol,:);
% met_y09 = objectives(bestSol,:);

[params_y09, met_y09] = sce_ua2(obj, x0,ngen, npop, npar, mag, lb, ub);

disp(strcat('Elapsed time: ', string(toc), 's'))

[Y_y09]= yates09(Setup, exp(params_y09));

% %% ShoreFor (Davidson et al. 2013)
% 

npar = 3;

x0 = log([100 3.7057e-05 64]);
lb = log([60 0.5e-5 50]);
ub = log([200 5e-3 80]);

disp("Starting ShoreFor...")
obj = @(X) Objective("SF", Setup, X, "RMSE", ENS);
% obj = @(X) MultiObjective("MD", Setup, X, ENS);
tic;

% [population, objectives] = nsga2(obj, npar, ngen, npop, ub, lb, nobj);
% 
% [~, bestSol] = min(sum(objectives,2));
% 
% params_sf = population(bestSol,:);
% met_sf = objectives(bestSol,:);

[params_sf, met_sf] = sce_ua2(obj, x0,ngen, npop, npar, mag, lb, ub);

disp(strcat('Elapsed time: ', string(toc), 's'))

[Y_sf]= shorefor(Setup, exp(params_sf));