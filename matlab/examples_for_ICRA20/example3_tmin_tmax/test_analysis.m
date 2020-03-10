clc
clear all

%% Grid
grid_min = [-2; -3; -0.1];  % Lower corner of computation domain
grid_max = [2; 3; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Create human dynamical system
% Case 1: Human's beta* = 0 (optimal)
%         Robot's Prior(beta = 0) = 0.1 (robot thinks human is random)
%          
%         Analysis 1.1: min_u computes the min time it takes robot to
%                       realize human is actually optimal
%         Analysis 1.2: max_u computes the max time it takes robot to
%                       realize human is actually optimal
%
% Case 2: Human's beta* = 1 (random)
%         Robot's Prior(beta = 0) = 0.9 (robot thinks human is optimal)
%
%         Analysis 2.1: min_u computes the min time it takes robot to
%                       realize human is actually random
%         Analysis 2.2: max_u computes the max time it takes robot to
%                       realize human is actually random

% The dynamics are too stiff for the HJI reachability for the computation of 
% Tmax so we use a smaller grid and a lower velocity to increase the accuracy.  
% Velocity
v = 0.1;

% Control bounds
uRange = [-pi+1e-2; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 1;

% Number of discrete controls
numCtrls = 11;

% Variance in normal distribution
sigma = pi/8;

% Threshold to determine likely controls
uThresh = 0.1;  

% Are we using dynamic of static beta model?
betaModel = 'static';

% Known human goal locations. 
goals = {[2,2], [2,-2]}; %{[1, tan(pi/4)], [1, -tan(pi/4)]};

% We have no dynamic beta parameters
extraArgs = [];

% Tolerance for how "sufficiently high" the probability needs to be
tol = 0.1;

% ---- Setup for Case 1 ---- %
% start with high prior on beta=1, but true is beta = 0 (gaussian)
Pgoal1 = 0.1; 
trueGoal = 1;
centerGoal = 1;

% % For Analysis 1.1 (min time)
% uMode = 'min';

% For Analysis 1.2 (max time)
uMode = 'max';
% -------------------------- %

% % ---- Setup for Case 2 ---- %
% % start with high prior on beta=0, but true is beta = 1 (rand)
% Pgoal1 = 0.9; 
% trueGoal = 2;
% centerGoal = 0;
% 
% % For Analysis 2.1 (min time)
% uMode = 'min';
% 
% % For Analysis 2.2 (max time)
% uMode = 'max';  
% % -------------------------- %

% Setup dynamical system
x0 = [0.; 0.; Pgoal1];
human = Goal1or0Human(x0, v, trueGoal, uRange, gamma, goals, sigma, ...
    uThresh, numCtrls, betaModel, extraArgs);

%% Setup target set
% Target set is centered at the true beta value
xyoffset = 0.1;
poffset = 0.01;
center = [0; 0; centerGoal];
widths = [(grid_max(1) - grid_min(1)) - xyoffset; ...
          (grid_max(2) - grid_min(2)) - xyoffset; 
          tol - poffset];
% widths = [2; 2; tol - poffset];
data0 = shapeRectangleByCenter(g, center, widths);

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
schemeData.hamFunc = @gaussianTwoGoalHuman_ham;
schemeData.partialFunc = @gaussianTwoGoalHuman_partial;

%% Setup value function computation params.
% HJIextraArgs.visualize = true; %show plot
HJIextraArgs.quiet = true;
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [-4 4 -4 4 -0.1 1.1];
HJIextraArgs.visualize.xTitle = '$p^x$';
HJIextraArgs.visualize.yTitle = '$p^y$';
HJIextraArgs.visualize.zTitle = '$P(Goal = 1)$';
HJIextraArgs.visualize.fontSize = 15;
HJIextraArgs.stopInit = x0;
%HJIextraArgs.visualize.camlightPosition = [0 0 0];

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid
HJIextraArgs.ignoreBoundary = 0; 

% Set the l(x) or target set in the min with l(x) formulation. 
%HJIextraArgs.targets = data0;

%uncomment if you want to see a 2D slice
HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
HJIextraArgs.visualize.plotData.projpt = 0.1; %project pt
HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%minWith = 'set';
minWith = 'zero';
%minWith = 'minVwithL';

%% Debugging
fprintf('------ Prediction Analysis Setup -------\n');
fprintf('   true human goal: %d\n', trueGoal);
fprintf('   [x, y, P(goal=1)]: [%d, %d, %d]\n', x0(1), x0(2), x0(3));
fprintf('   uMode: %s\n', uMode);
fprintf('--------------------------------------\n');


%% Solve it!
[data, tau2, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);


% Compute the optimal control
deriv = computeGradients(g, data(:, :, :, end));
%  Compute the optimal over all states and discretized controls.
p1 = repmat(deriv{1}, [1, 1, 1, human.numCtrls]);
p2 = repmat(deriv{2}, [1, 1, 1, human.numCtrls]);
p3 = repmat(deriv{3}, [1, 1, 1, human.numCtrls]);
pdot_f = p1.*human.xdot{1} + p2.*human.xdot{2} + p3.*human.xdot{3};
% Compute the optimizing control
if strcmp(uMode, 'min')
	[~, indices] = min(pdot_f, [], 4);
else
	[~, indices] = max(pdot_f, [], 4);
end
% Compute the optimal control
optCtrl_index = indices(41, 41, 10); %Optimal control index at x0
optCtrl = human.likelyCtrls{optCtrl_index}(41, 41, 10); % The actula control
probCtrl = normpdf((optCtrl-pi/4)/sigma);

%% Grab the min/max time. 
% if trueBeta == 0
%     % If the true beta = 0, we want P(beta=0) >= (1-tol)
%     lowIdx = PbetaToGrid(g, 1-tol);
%     upIdx = PbetaToGrid(g, 1);
% else
%     % If the true beta = 1, we want P(beta=0) <= tol
%    lowIdx = PbetaToGrid(g, 0);
%    upIdx = PbetaToGrid(g, tol);
% end

% tmin = -Inf;
% for t=1:length(tau2)
%     v = eval_u(g, data(:,:,:,t), x0); 
%     if v <= 0
%         tmin = (t-1)*dt;
%         break;
%     end
% end
% 
% fprintf("Minimum time it takes to realize beta=%d is %f\n", trueBeta, tmin);
% 
% %% Converts from (pbeta) state to (i) grid index.
% function PbetaIdx = PbetaToGrid(grid, pb)
%     error = abs(grid.xs{3} - pb);
%     [~,idx] = min(error(:));
%     [~,~,PbetaIdx] = ind2sub(size(error),idx);
% end