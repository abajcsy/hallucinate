function params = frsTwoGoalPrecompute()
%% Predictor: Create human dynamical system for prediction
vH = 0.6;                           % Velocity
uHRange = [-pi+1e-2; pi];           % Control bounds
gamma = 1;                          % gamma in continuous-time P(beta = 0) dynamics

numCtrls = 11;                      % Number of discrete controls
uHThresh = 0.05;                    % Threshold to determine likely controls

sigma = pi/4;                               % Variance in normal distribution
goals = {[2,2], [2,-2]};  % Known human goal locations. 

betaModel = 'static';               % Are we using dynamic of static beta model?
extraArgs = [];                     % (No) Dynamic beta parameters

Pgoal1 = 0.5;                       % Prior over goal=goal1

params.xH0 = [0; 0];                % Initial physical state of human    
params.z0 = [params.xH0; Pgoal1];   % Initial condition for reachability

% Create human prediction model.
params.humanType = "twoGoal";
params.humanModel = GaussianTwoGoalHuman(params.z0, vH, uHRange, gamma, ...
    goals, sigma, uHThresh, numCtrls, betaModel, extraArgs);

% Setup custom hamiltonians.
params.hamFunc = @gaussianTwoGoalHuman_ham;
params.partialFunc = @gaussianTwoGoalHuman_partial;

%% Environment Params (meters).
params.lowEnv = [-4;-4];
params.upEnv = [4;4];

%% Forward Reachability: Prediction Computation Params.

% Discretization of x,y, and P(beta=0) space.
predGridLow = [params.lowEnv; -0.1];
predGridUp = [params.upEnv; 1.1];
params.N = [81; 81; 41];
params.predGrid = createGrid(predGridLow, predGridUp, params.N);

% Let the human have access to the grid for debugging.
params.humanModel.setGrid(params.predGrid);

% Pre-compute the optimal control over the entire state-space.
params.humanModel.computeUOptGoals(params.predGrid.xs);

% Pre-compute the likely controls and dynamics over the entire state-space.
params.humanModel.computeUAndXDot(params.predGrid.xs);

% Prediction time horizon and discretization.
params.tMin = 0;
params.tMax = 6.2;
params.dt = 0.05;

% Target set radius.
params.targetRad = 0.1;

% Setup for the control 
params.uMode = 'max'; 
params.minWith = 'set'; % minwith = 'zero' gives us tube.
params.quiet = true; % runs in quiet computation mode, for efficiency.

end