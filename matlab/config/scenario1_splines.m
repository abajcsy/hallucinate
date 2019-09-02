function params = scenario1_splines()
%% Scenario 1
%  Robot needs to cross paths with human to get to its goal. If it trusts
%  that the human will continue to move forward, it will overconfidently
%  plan to cut in front of the human and end up in a collision. 

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

params.xH0 = [-0.5; -0.5];            % Initial physical state of human    
params.z0 = [params.xH0; Pbeta0];   % Initial condition for reachability

% Create human prediction model.
params.humanType = "gaussian";
params.humanModel = GaussianHuman(params.z0, vH, uHRange, gamma, K, m, ...
                             sigma, uHThresh, numCtrls, betaModel, extraArgs);

% Setup custom hamiltonians.
params.hamFunc = @gaussianHuman_ham;
params.partialFunc = @gaussianHuman_partial;

% Path to FRS directory (assuming the root directory is set to the project
% root directory in MATLAB).
params.pathToFRSDir = './matlab/frs';

%% Human Simluator: Parameters for simulating human measurements. 

x0Sim = params.xH0;
vHSim = vH;
ctrls = {0, pi/4, pi/2, (3*pi)/4, pi};
times = {[0,0.5], [0.5,0.7], [0.7,1.0], [1.0,1.3], [1.3,4]};

params.simHuman = FixedTrajHuman(x0Sim, vHSim, ctrls, times);

%% Environment Params (meters).
params.lowEnv = [-2;-2];
params.upEnv = [2;2];

rect1 = {params.lowEnv, [params.upEnv(1); params.lowEnv(2)+0.5]};
rect2 = {[params.lowEnv(1); params.upEnv(2)-0.5], params.upEnv};
params.staticObsBounds = {rect1, rect2}; % has lower and upper bounds for 2 static obstacles.

%% Forward Reachability: Prediction Computation Params.

% Discretization of x,y, and P(beta=0) space.
predGridLow = [params.lowEnv; -0.1];
predGridUp = [params.upEnv; 1.1];
params.N = [81; 81; 41];
params.predGrid = createGrid(predGridLow, predGridUp, params.N);

% Pre-compute the likely controls and dynamics over the entire state-space.
params.humanModel.computeUAndXDot(params.predGrid.xs);

% Prediction time horizon and discretization.
params.tMin = 0;
params.tMax = 2;
params.dt = 0.05;

% Target set radius.
params.targetRad = 0.1;

% Setup for the control 
params.uMode = 'max'; 
params.minWith = 'set'; % minwith = 'zero' gives us tube.
params.quiet = true; % runs in quiet computation mode, for efficiency.

% Color of predictions for plotting.
params.predColor = [99., 180., 255.]/255.;

%% Robot: Planning Params.

% Discretization in (x,y,theta,v,time)
params.xDisc = 0.15;
params.yDisc = 0.15;
% params.xDisc = 0.25;
% params.yDisc = 0.25;
params.thetaDisc = pi / 45;
% params.vDisc = 0.5;
params.vDisc = 0.1;
params.tDisc = 0.5;

% Heuristic and goal parameters.
% params.heurWeight = 1.25;
% params.heurWeight = 1;
params.heurWeight = 2;
params.goalTol = 0.3;

params.minEdgeTimeSteps = 2;
params.maxEdgeTimeSteps = 2;

% Setup the state bounds.
params.xBounds = [params.lowEnv(1), params.upEnv(1)];
params.yBounds = [params.lowEnv(2), params.upEnv(2)];
params.vBounds = [0, 0.6];
params.thetaBounds = [-2*pi, 2*pi];
params.timeBounds = [0, 15]; % Planning horizon.

% Setup initial state of planner and goal state of planner.
params.xR0 = [0.1; -1; (3*pi)/4; 0.1];
params.goalRXY = [-1; 1];

%% Robot: Dynamical System Params.
wMax = 1;
%aRange = [0,0.1];
vRange = params.vBounds;
dMax = 0;
params.simRobot = Plane(params.xR0(1:3), wMax, vRange, dMax);
% params.simRobot = DubinsCar(params.xR0(1:3), wMax, speed, dMax);
% params.simRobot = Unicycle4DRobot(params.xR0(1:4), wMax, aRange, vRange);

%% Simulation Params.
% Timestep for computation and simulation.
params.simDt = 0.05;
% params.simDt = 0.045;
% params.simDt = 0.1;
% params.simDt = 1;
% params.simDt = 0.5;
% params.simDt = 0.11;
params.T = 200;

% Number of steps after which to replan.
% params.replanAfterSteps = 1;
% params.replanAfterSteps = 3;
params.replanAfterSteps = params.T + 1;

% If set to true, uses the control from the trajectory; otherwise, uses the
% state directly.
params.trajUseControl = false;

end