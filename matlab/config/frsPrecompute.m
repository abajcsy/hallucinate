function params = frsPrecompute()
%% Predictor: Create human dynamical system for prediction
vH = 0.6;                           % Velocity
uHRange = [-pi+1e-2; pi];            % Control bounds
gamma = 1;                          % gamma in continuous-time P(beta = 0) dynamics

K = [0, 0];                         % Control gains
m = 0;

numCtrls = 11;                      % Number of discrete controls
sigma = 0.1;                        % Variance in normal distribution
uHThresh = 0.05;                    % Threshold to determine likely controls
Pbeta0 = 0.5;                       % Prior over beta = 0 (rational human)

betaModel = 'static';               % Are we using dynamic of static beta model?
extraArgs.alpha = 0.5;              % (dynamic beta) parameter
extraArgs.DeltaB0 = 0.5;            % (dynamic beta) mixing distribution 

params.xH0 = [0; 0];                % Initial physical state of human    
params.z0 = [params.xH0; Pbeta0];   % Initial condition for reachability

% Create human prediction model.
params.humanType = "gaussian";
params.humanModel = GaussianHuman(params.z0, vH, uHRange, gamma, K, m, ...
                             sigma, uHThresh, numCtrls, betaModel, extraArgs);

% Setup custom hamiltonians.
params.hamFunc = @gaussianHuman_ham;
params.partialFunc = @gaussianHuman_partial;

%% Environment Params (meters).
params.lowEnv = [-4;-4];
params.upEnv = [4;4];

rect1 = {params.lowEnv, [params.upEnv(1); params.lowEnv(2)+0.5]};
rect2 = {[params.lowEnv(1); params.upEnv(2)-0.5], params.upEnv};
params.staticObsBounds = {rect1, rect2}; % has lower and upper bounds for 2 static obstacles.

%% Forward Reachability: Prediction Computation Params.

% Discretization of x,y, and P(beta=0) space.
predGridLow = [params.lowEnv; -0.1];
predGridUp = [params.upEnv; 1.1];
params.N = [121; 121; 41];
params.predGrid = createGrid(predGridLow, predGridUp, params.N);

% Pre-compute the likely controls and dynamics over the entire state-space.
params.humanModel.computeUAndXDot(params.predGrid.xs);

% Prediction time horizon and discretization.
params.tMin = 0;
params.tMax = 6;
params.dt = 0.05;

% Target set radius.
params.targetRad = 0.2;

% Setup for the control 
params.uMode = 'max'; 
params.minWith = 'set'; % minwith = 'zero' gives us tube.
params.quiet = true; % runs in quiet computation mode, for efficiency.

end