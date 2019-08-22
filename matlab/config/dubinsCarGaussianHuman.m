function params = dubinsCarGaussianHuman()

%% Create human dynamical system for *prediction*
v = 0.6;                        % Velocity
uRange = [-pi; pi];             % Control bounds
gamma = 1;                      % gamma in continuous-time P(beta = 0) dynamics

K = [0, 0];                     % Control gains
m = 0;

numCtrls = 10;                  % Number of discrete controls
sigma = 0.1;                    % Variance in normal distribution
uThresh = 0.05;                 % Threshold to determine likely controls
DeltaB0 = 0.5;                  % Distribution in HMM
Pbeta0 = 0.001;                   % Setup dynamical system

params.x0 = [0; 0; Pbeta0];     % Initial condition

% Create human model.
params.humanType = "gaussian";
params.humanModel = GaussianHuman(params.x0, v, uRange, gamma, K, m, ...
                             sigma, uThresh, DeltaB0, numCtrls);

% Setup custom hamiltonians.
params.hamFunc = @gaussianHuman_ham;
params.partialFunc = @gaussianHuman_partial;

%% Simulated human params.

% beta = 1 --> irrational
% beta = 0 --> rational
trueBeta = 1; 
mu = K*params.x0(1:2) + m;
params.simHuman = SimFixedBetaGaussianHuman(params.x0(1:2), v, ...
                                            mu, sigma, uRange, trueBeta);

%% Environment Params (meters).
params.lowEnv = [-2;-2];
params.upEnv = [2;2];

%% Prediction Computation Params.

% Discretization of x,y, and P(beta=0) space.
gridLow = [params.lowEnv; -0.1];
gridUp = [params.upEnv; 1.1];
params.N = [81; 81; 41];
params.predGrid = createGrid(gridLow, gridUp, params.N);

% Prediction time horizon and discretization.
params.tMin = 0;
params.tMax = 2;
params.dt = 0.05;

% Target set radius.
params.targetR = 0.2;

% Setup for the control 
params.uMode = 'max'; 
params.minWith = 'set'; % minwith = 'zero' gives us tube.
params.quiet = true; % runs in quiet computation mode, for efficiency.

%% Planning Params.
params.xDisc = 0.15;
params.yDisc = 0.15;
params.thetaDisc = pi / 45;
params.vDisc = 0.5;
params.tDisc = 0.25;
params.heurWeight = 1.25;
params.goalTol = 0.2;

%% Robot Dynamical System Params.
params.robotx0 = [0; 0; 0; 1; 0];
params.robotGoalXY = [1.8; 1.8];

%% Simulation Params.
% Timestep for computation and simulation.
params.simDt = 0.05;
params.T = 2000; 

end