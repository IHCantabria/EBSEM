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


%% Apply models

params_md = [4.3685,   -9.5850,   -7.0647];
params_y09 = [-2.1862,    2.2438,   -5.6121,   -5.6156];
params_sf = [4.6040,  -10.2061,    4.1673];

Y_md = millerDean04(Setup, exp(params_md));
Y_y09 = yates09(Setup, exp(params_y09));
Y_sf = shorefor(Setup, exp(params_sf));

figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.2]);
scatter(ENS.time, ENS.Yobs, 1, 'MarkerEdgeColor', [0.5, 0.5, 0.5], 'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on;
plot(time, Y_md, 'b', 'LineWidth', 0.5);
plot(time, Y_y09, 'r', 'LineWidth', 0.5);
plot(time, Y_sf, 'k', 'LineWidth', 0.5);
hold off;
ylabel('Y [m]', 'FontSize', 12);
legend('Observed data', 'Miller and Dean 2004', 'Yates et al. 2009', 'Davidson et al. 2013', 'FontSize', 8, 'Location', 'north', 'NumColumns', 4);
xlim([time(1), time(end)]);
grid on;
set(gca, 'GridLineStyle', '--', 'LineWidth', 0.5);
datetick
