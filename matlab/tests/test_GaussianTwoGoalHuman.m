clc
clear all

%% Create human dynamical system

% Velocity
v = 0.6;

% Control bounds
uRange = [-pi+1e-2; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 1;

% Number of discrete controls
numCtrls = 11;

% Threshold to determine likely controls
uThresh = 0.05; 

% Variance in normal distributions
sigma = pi/4;

% Known human goal locations. 
goals = {[1, tan(pi/4)], [1, -tan(pi/4)]};

% Are we using dynamic of static beta model?
betaModel = 'static';

% (No) Dynamic beta parameters
extraArgs = [];

% Setup dynamical system
Pgoal1 = 0.5; 
x0 = [0; 0; Pgoal1];
human = GaussianTwoGoalHuman(x0, v, uRange, gamma, goals, sigma, uThresh, numCtrls, ...
    betaModel, extraArgs);

%% Grid
grid_min = [-2; -2; -0.1];  % Lower corner of computation domain
grid_max = [2; 2; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Let the human have access to the grid for debugging.
human.setGrid(g);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUOptGoals(g.xs);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% target set
% Findings: (1) Increasing R reduces the gap between the two extremes so it 
% could just be numerical issues. (2) Increasing the grid resolution for
% R=0.1 also reduces the gap between the two extremes so that is other
% evidence that it is just numerical error.
R = 0.1;
data0 = shapeSphere(g, x0, R);

%% time vector
t0 = 0;
tMax = 3;
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
schemeData.hamFunc = @gaussianTwoGoalHuman_ham;
schemeData.partialFunc = @gaussianTwoGoalHuman_partial;

%% Compute value function
% HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [-2 2 -2 2 -0.1 1.1];
HJIextraArgs.visualize.xTitle = '$p^x$';
HJIextraArgs.visualize.yTitle = '$p^y$';
HJIextraArgs.visualize.zTitle = '$P(\beta = 0)$';
HJIextraArgs.visualize.fontSize = 15;
%HJIextraArgs.visualize.camlightPosition = [0 0 0];

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid
HJIextraArgs.ignoreBoundary = 0; 

%uncomment if you want to see a 2D slice
HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
HJIextraArgs.visualize.plotData.projpt = {'min'}; %project pt
HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%HJIextraArgs.targets = data0;

minWith = 'set';
%minWith = 'zero';
%minWith = 'minVwithL';
[data, tau2, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);

