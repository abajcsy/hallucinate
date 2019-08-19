clc
clear all

%% NOTE:
% In HJIPDESolve, uncomment the following lines:
%   if isfield(schemeData, 'dynSys')
%     schemeData.hamFunc = @genericHam;
%     schemeData.partialFunc = @genericPartial;
%   end

%% Create human dynamical system

% Velocity
v = 1;

% Control bounds
uRange = [-pi; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 0;

% Control gains
K = [0, 0];
m = 0;

% Number of discrete controls
numCtrls = 10;

% Variance in normal distribution
sigma = 0.1;

% Threshold to determine likely controls
uThresh = 0.05; 

% Distribution in HMM
DeltaB0 = 0.5; 

% Setup dynamical system
Pbeta0 = 0.9; 
x0 = [0; 0; Pbeta0];
human = GaussianHuman(x0, v, uRange, gamma, K, m, sigma, uThresh, DeltaB0, numCtrls);

%% Grid
grid_min = [-2; -2; -0.1]; % Lower corner of computation domain
grid_max = [2; 2; 1.1];    % Upper corner of computation domain
N = [81; 81; 41];         % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% target set
R = 0.2;
data0 = shapeSphere(g, x0, R);

%% time vector
t0 = 0;
tMax = 2;
dt = 0.05;
tau = t0:dt:tMax;
uMode = 'max';

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'forward';
schemeData.hamFunc = @gaussianHuman_ham;
schemeData.partialFunc = @gaussianHuman_partial;

%% Compute value function
%HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = false;
HJIextraArgs.visualize.viewAxis = [-2 2 -2 2 0 1.1];
%az = -10;
%el = 33;
%HJIextraArgs.visualize.viewAngle = [az, el];
HJIextraArgs.visualize.xTitle = "$p^x$";
HJIextraArgs.visualize.yTitle = "$p^y$";
HJIextraArgs.visualize.zTitle = "$P(\beta = 0)$";
HJIextraArgs.visualize.fontSize = 15;

%uncomment if you want to see a 2D slice
HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
HJIextraArgs.visualize.plotData.projpt = {'min'}; %project pt
HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

minWith = 'set';
%minWith = 'zero';
[data, tau2, ~] = ...
  HJIPDE_solve(data0, tau, schemeData, minWith, HJIextraArgs);

