clc
clear all

%% Create human dynamical system
v = 1;
uRange = [-pi; pi];
alpha = 1;
K = [0, 0];
m = 0;
uThresh = 0.05;
betaPrior = 0.5; % Delta_0 = P^-_0(\beta = 0)
x0 = [0; 0; betaPrior];
human = FeedbackHuman(x0, v, uRange, alpha, K, m, uThresh, betaPrior);

%% Grid
grid_min = [-2; -2; -0.1]; % Lower corner of computation domain
grid_max = [2; 2; 1.1];    % Upper corner of computation domain
N = [81; 81; 41];         % Number of grid points per dimension
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
schemeData.tMode = 'forward';
schemeData.hamFunc = @feedbackHuman_ham;
schemeData.partialFunc = @feedbackHuman_partial;

%% Compute value function
%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = false;
HJIextraArgs.visualize.viewAxis = [-2 2 -2 2 0.1 1.1];
HJIextraArgs.visualize.xTitle = "$x$";
HJIextraArgs.visualize.yTitle = "$y$";
HJIextraArgs.visualize.zTitle = "$P(\beta = 0)$";
HJIextraArgs.visualize.fontSize = 15;

%HJIextraArgs.makeVideo = 1;
%HJIextraArgs.videoFilename = "frs_beta.mp4";

%uncomment if you want to see a 2D slice
% HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
% HJIextraArgs.visualize.plotData.projpt = [0]; %project pt
% HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

minWith = 'zero';

%[data, tau, extraOuts] = ...
% HJIPDE_solve(data0, tau, schemeData, minWith, extraArgs)
[data, tau2, ~] = ...
  HJIPDE_solve(data0, tau, schemeData, minWith, HJIextraArgs);

