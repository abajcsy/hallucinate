clc
clear all

%% Create human dynamical system

% Velocity
v = 0.6;

% Control bounds
uRange = [-pi+1e-2; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 0.01;

% Control gains
K = [0, 0];
m = 0;

% Number of discrete controls
numCtrls = 20;

delta_t = 1;

% Threshold to determine likely controls
uThresh = 0.5;

% Are we using dynamic of static beta model?
betaModel = 'static';

theta = [0 1];

% Dynamic beta parameters
extraArgs.alpha = 0.5;
extraArgs.DeltaB0 = 0.5; 

% Setup dynamical system
Pbeta0 = 0.8; 
x0 = [0; 0; Pbeta0];
human = BoltzmannHuman(x0, v, uRange, gamma, K, m, theta, delta_t, uThresh, ...
    numCtrls, betaModel, extraArgs);

%% Grid
grid_min = [-2; -2; -0.1];  % Lower corner of computation domain
grid_max = [2; 2; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Compute P(u_t|x_t,\beta=0)
u = 0; % Control
x = g.xs; 

lb = human.uRange(1);
ub = human.uRange(2);
inc = (ub-lb)/(human.numCtrls-1);

intControls = 0.0;
for i = 1:human.numCtrls

    % Find control
    u = lb + (i-1)*inc;

    % Find next x by forward euler
    x1 = x{1} + human.deltaT .* human.v .* cos(u);
    x2 = x{2} + human.deltaT .* human.v .* sin(u);

    % Evaluate distance of next x to goal theta under L2 norm
    nrm = ((x1 - human.theta(1)).^2 + (x2 - human.theta(2)).^2).^(0.5);

    % Calculate value in summation: e^{-\| (x_t + \Delat t f(x_t,u_t)) - \theta \|_2}
    val = exp(1).^(-1.*nrm);

    % Add to running value of summation
    intControls = intControls + inc*val;
end

% x_{t+1} = x_t + \Delta t * f(x_t)
x1 = x{1} + human.deltaT .* human.v .* cos(u);
x2 = x{2} + human.deltaT .* human.v .* sin(u);

nrm = ((x1 - human.theta(1)).^2 + (x2 - human.theta(2)).^2).^(0.5);
val = exp(1).^(-1.*nrm);

prob_u = val ./ intControls;
prob_u_high = prob_u(:,:,38); % high posterior prob of ~1.0

heatmap(prob_u_high)
