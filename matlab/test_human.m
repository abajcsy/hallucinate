clc
clear all

%% Create human dynamical system
v = 1;
x0 = [0; 0; 0.5];
uRange = [-pi; pi];
alpha = 0.02;
K = [0, 0];
m = 0;
uThresh = 0.1;
betaPrior = 0.99;
human = FeedbackHuman(x0, v, uRange, alpha, K, m, uThresh, betaPrior);

%% Grid
grid_min = [-1; -1; -0.1]; % Lower corner of computation domain
grid_max = [1; 1; 1.1];    % Upper corner of computation domain
N = [41; 41; 41];         % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% target set
R = 0.1;
% data0 = shapeCylinder(grid,ignoreDims,center,radius)
data0 = shapeSphere(g, x0, R);

%% time vector
t0 = 0;
tMax = 5;
dt = 0.05;
tau = t0:dt:tMax;
uMode = 'min';

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'high'; %set accuracy
schemeData.uMode = uMode;

%% Compute value function
%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update

% uncomment if you want to see a 2D slice
%HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
%HJIextraArgs.visualize.plotData.projpt = [0]; %project at theta = 0
%HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

minWith = 'zero';

%[data, tau, extraOuts] = ...
% HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
[data, tau2, ~] = ...
  HJIPDE_solve(data0, tau, schemeData, minWith, HJIextraArgs);

