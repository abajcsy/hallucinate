clc
clf
clear all

%% Grid
grid_min = [-5; -5; -0.1];  % Lower corner of computation domain
grid_max = [5; 5; 1.1];     % Upper corner of computation domain
N = [41; 41; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Create human dynamical system
% Case 1: Human's lambda* = 1 (goal 1)
%         Robot's Prior(lambda = 1) = 0.1 (robot thinks human is going to goal 1)
%          
%         Analysis 1.1: min_u computes the min time it takes robot to
%                       realize human is actually going to goal 1
%         Analysis 1.2: max_u computes the max time it takes robot to
%                       realize human is actually going to goal 1
%
% Case 2: Human's lambda* = 2 (goal 2)
%         Robot's Prior(lambda = 0) = 0.9 (robot thinks human is going to goal 2)
%
%         Analysis 2.1: min_u computes the min time it takes robot to
%                       realize human is actually going to goal 1
%         Analysis 2.2: max_u computes the max time it takes robot to
%                       realize human is actually going to goal 1

% Velocity
v = 0.2;

% Control bounds
uRange = [-pi+1e-2; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 1;

% Number of discrete controls
numCtrls = 31;

% Variance in normal distributions
sigma = pi/8;

% Known human goal locations. 
goals = {[2,-2], [2,2]};
%goals = {[2,-0.5], [2,0.5]};
%goals = {[2, 0], [3, 0]};

% Threshold to determine likely controls
uThresh = 0.1;  

% Are we using dynamic of static parameter model?
betaModel = 'static';

% We have no dynamic beta parameters
extraArgs = [];

% Tolerance for how "sufficiently high" the probability needs to be for us
% to be confident in the model. 
tol = 0.1;

% ---- Setup for Case 1 ---- %
% start with high prior on beta=1, but true is goal = 1 
Pgoal1 = 0.5; 
trueGoalIdx = 1;
centerPgoal1 = 1;

% For Analysis 1.1 (min time)
uMode = 'min';

% Setup dynamical system
x0 = [0; 0; Pgoal1];
human = GaussianG1orG2Human(x0, v, trueGoalIdx, uRange, gamma, goals, ...
                        sigma, uThresh, numCtrls, betaModel, extraArgs);

%% Setup target set
% Target set is centered at the true beta value
xyoffset = 0.1;
poffset = 0.01;
center = [0; 0; centerPgoal1];
widths = [(grid_max(1) - grid_min(1)) - xyoffset; ...
          (grid_max(2) - grid_min(2)) - xyoffset; 
          tol - poffset];
data0 = shapeRectangleByCenter(g, center, widths);

%% Let the human have access to the grid for debugging.
human.setGrid(g);

%% Pre-compute the optimal control over the entire state-space.
human.computeUOptGoals(g.xs);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% time vector
t0 = 0;
tMax = 100;
dt = 0.1;
tau = t0:dt:tMax;

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'backward';
schemeData.hamFunc = @gaussianG1orG2Human_ham;
schemeData.partialFunc = @gaussianG1orG2Human_partial;

%% Setup value function computation params.
HJIextraArgs.quiet = false;
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [grid_min(1) grid_max(1) grid_min(2) grid_max(2) -0.1 1.1];
HJIextraArgs.visualize.xTitle = '$x$';
HJIextraArgs.visualize.yTitle = '$y$';
HJIextraArgs.visualize.zTitle = '$P(goal = 1)$';
HJIextraArgs.visualize.fontSize = 15;
HJIextraArgs.stopInit = x0;

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid
HJIextraArgs.ignoreBoundary = 0; 

% Min with zero is more robust to "vanishing volume" issue.
%minWith = 'set';
minWith = 'zero';

%% Debugging
fprintf('------ Prediction Analysis Setup -------\n');
fprintf('   true human goal: %d\n', trueGoalIdx);
fprintf('   [x, y, P(goal=1)]: [%d, %d, %d]\n', x0(1), x0(2), x0(3));
fprintf('   uMode: %s\n', uMode);
fprintf('--------------------------------------\n');


%% Solve it!
[valueFuns, times, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);

%% Grab the min time. 
tmin = times(end);
tminIdx = length(times);

% tmin = -Inf;
% tminIdx = 0;
% for t=1:length(times)
%     v = eval_u(g, valueFuns(:,:,:,t), x0); 
%     if v <= 0
%         tmin = (t-1)*dt;
%         tminIdx = t;
%         break;
%     end
% end

fprintf('Minimum time it takes to realize goal=%d is %f\n', trueGoalIdx, tmin);

%% Get the sequence of optimal controls.
uopt = GetOptControls(x0, g, valueFuns, times(1:tminIdx), human, uMode);

%% Simulate and visualize optimal control and states.
plotStateTraj(uopt, human, times, tminIdx, goals, trueGoalIdx, grid_min, grid_max);