clc
clear all

%% Create human dynamical system

% Velocity
v = 0.2;

% Control bounds
uRange = [-pi+1e-2; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 1;

% Control gains
K = [0, 0];
m = 0;

% Number of discrete controls
numCtrls = 21;

% Variance in normal distribution
sigma = 0.1;

% Threshold to determine likely controls
uThresh = 0.05; 

% Are we using dynamic of static beta model?
betaModel = 'static';

% We have no dynamic beta parameters
extraArgs = [];

% Tolerance for how "sufficiently high" the probability needs to be
tol = 0.1;

% Setup dynamical system
% Some options:
%   Pbeta0 = 0.9, trueBeta = 1
%       start with high prior on beta=0, but true is beta = 1 (rand)
%   Pbeta = 0.1, trueBeta = 0
%       start with high prior on beta=1, but true is beta = 0 (gaussian)
Pbeta0 = 0.1; 
x0 = [0; 0; Pbeta0];
trueBeta = 0;
human = Beta1or0Human(x0, v, trueBeta, uRange, gamma, K, m, sigma, ...
    uThresh, numCtrls, betaModel, extraArgs);

% Choose 'max' if computing MIN Delta T
% Choose 'min' if computing MAX Delta T
uMode = 'max';

%% Grid
grid_min = [-2; -2; -0.1];  % Lower corner of computation domain
grid_max = [2; 2; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% target set
% Findings: (1) Increasing R reduces the gap between the two extremes so it 
% could just be numerical issues. (2) Increasing the grid resolution for
% R=0.1 also reduces the gap between the two extremes so that is other
% evidence that it is just numerical error.
R = 0.2;
data0 = shapeSphere(g, x0, R);

%% time vector
t0 = 0;
tMax = 5;
dt = 0.1;
tau = t0:dt:tMax;

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'forward';
schemeData.hamFunc = @beta1or0Human_ham;
schemeData.partialFunc = @beta1or0Human_partial;

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

% %uncomment if you want to see a 2D slice
% HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
% HJIextraArgs.visualize.plotData.projpt = {'min'}; %project pt
% HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%HJIextraArgs.targets = data0;

%minWith = 'set';
minWith = 'zero';
%minWith = 'minVwithL';
[data, tau2, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);

if trueBeta == 0
    % If the true beta = 0, we want P(beta=0) >= (1-tol)
    lowIdx = PbetaToGrid(g, 1-tol);
    upIdx = PbetaToGrid(g, 1);
else
    % If the true beta = 1, we want P(beta=0) <= tol
   lowIdx = PbetaToGrid(g, 0);
   upIdx = PbetaToGrid(g, tol);
end

tmin = -Inf;
for t=1:length(tau2)
    v = data(:,:,:,t);
    B = (v(:,:,lowIdx:upIdx) <= 0);
    total = sum(B(:));
    if total > 0
        tmin = t*dt;
        break;
    end
end

fprintf("Minimum time it takes to realize beta=%d is %f\n", trueBeta, tmin);

%% Converts from (pbeta) state to (i) grid index.
function PbetaIdx = PbetaToGrid(grid, pb)
    error = abs(grid.xs{3} - pb);
    [~,idx] = min(error(:));
    [~,~,PbetaIdx] = ind2sub(size(error),idx);
end