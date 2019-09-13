clc
clear all

%% Create human dynamical system

% Velocity
v = 0.6;

% Control bounds
uRange = [-pi+1e-2; pi];

% Control bounds for the robot
vRRange = [0., 0.6];
wRRange = [-1.1, 0.1];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 1;

% Control gains
K = [0, 0];
m = 0;

% Number of discrete controls
numCtrls = 11;

% Variance in normal distribution
sigma = 0.1;

% Threshold to determine likely controls
uThresh = 0.05; 

% Are we using dynamic of static beta model?
betaModel = 'static';

% Dynamic beta parameters
extraArgs.alpha = 0.5;
extraArgs.DeltaB0 = 0.5; 

% Setup dynamical system
Pbeta0 = 0.9; 
x0 = [0; 0; Pbeta0; 0];
human = GaussianHumanUnicycleRobot(x0, v, uRange, gamma, K, m, sigma, uThresh, numCtrls, ...
    betaModel, vRRange, wRRange, extraArgs);

%% Grid
grid_min = [-2; -2; -0.1; -pi];  % Lower corner of computation domain
grid_max = [2; 2; 1.1; pi];     % Upper corner of computation domain
N = [21; 21; 21; 21];           % Number of grid points per dimension
pdDims = 4; %Periodic dimension
g = createGrid(grid_min, grid_max, N, pdDims);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% target set
R = 0.3;
data0 = shapeCylinder(g, [3, 4], x0, R);

%% time vector
t0 = 0;
tMax = 2;
dt = 0.1;
tau = t0:dt:tMax;
uMode = 'max';
dMode = 'min';

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.dMode = dMode;
schemeData.tMode = 'backward';
schemeData.hamFunc = @gaussianHumanUnicycleRobot_ham;
schemeData.partialFunc = @gaussianHumanUnicycleRobot_partial;

%% Compute value function
% HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [-2 2 -2 2 -pi pi];
HJIextraArgs.visualize.xTitle = '$p^x$';
HJIextraArgs.visualize.yTitle = '$p^y$';
HJIextraArgs.visualize.zTitle = '$\theta$';
HJIextraArgs.visualize.fontSize = 15;
%HJIextraArgs.visualize.camlightPosition = [0 0 0];

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid
HJIextraArgs.ignoreBoundary = 0; 

%uncomment if you want to see a 2D slice
HJIextraArgs.visualize.plotData.plotDims = [1 1 0 1]; %plot x, y
HJIextraArgs.visualize.plotData.projpt = {0.1}; %project pt
HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%HJIextraArgs.targets = data0;

% minWith = 'set';
minWith = 'zero';
%minWith = 'minVwithL';
[data, tau2, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);

