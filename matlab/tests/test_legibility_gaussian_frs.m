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
numCtrls = 21;

% Variance in normal distributions
sigma = pi/4;

% Known human goal locations. 
goals = {[2,-2], [2,2]};
%goals = {[2,-1], [2,1]};
%goals = {[2, 0], [3, 0]};

% Threshold to determine likely controls
uThresh = 0.0;  

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
%uMode = 'min';
uMode = 'max';

% Setup dynamical system
x0 = [-2; 0; Pgoal1];
human = GaussianG1orG2Human(x0, v, trueGoalIdx, uRange, gamma, goals, ...
                        sigma, uThresh, numCtrls, betaModel, extraArgs);

%% Setup target set
% Target set is centered at the true beta value
% xyoffset = 0.1;
% poffset = 0.01;
% center = [0; 0; centerPgoal1];
% widths = [(grid_max(1) - grid_min(1)) - xyoffset; ...
%           (grid_max(2) - grid_min(2)) - xyoffset; 
%           tol - poffset];
% data0 = shapeRectangleByCenter(g, center, widths);
R = 0.2;
data0 = shapeSphere(g, x0, R);

%% Let the human have access to the grid for debugging.
human.setGrid(g);

%% Pre-compute the optimal control over the entire state-space.
human.computeUOptGoals(g.xs);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% time vector
t0 = 0;
tMax = 25;
dt = 0.1;
tau = t0:dt:tMax;

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'forward';
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

fprintf("Minimum time it takes to realize goal=%d is %f\n", trueGoalIdx, tmin);

%% Get the sequence of optimal controls.
uopt = GetOptControls(x0, g, valueFuns, times(1:tminIdx), human, uMode);

%% Simulate and visualize optimal control and states.
plotStateTraj(uopt, human, times, tminIdx, goals, trueGoalIdx, grid_min, grid_max);

%% Grab the sequence of optimal controls starting from a state x. 
%  Returns cell array of all the optimal controls.
function uopt = GetOptControls(x, grid, valueFuns, times, human, uMode)
    uopt = cell(1,length(times));
    for t=1:length(times)
        % Grab the derivative at all states.
        deriv = computeGradients(grid, valueFuns(:,:,:,t));

        % Value of the derivative at that particular state
        current_deriv = eval_u(grid, deriv, x);

        % Get the optimal control to apply at this state
        u = human.optCtrl(t, x, current_deriv, uMode); 
        uopt{t} = u;
    end
end

function plotStateTraj(uopt, human, times, tminIdx, goals, trueGoalIdx, grid_min, grid_max)
    figure(2);
    hold on

    % Keep track of all states so we can draw a line plot.
    X = zeros(1,length(times(1:tminIdx)));
    Y = zeros(1,length(times(1:tminIdx)));
    PG = zeros(1,length(times(1:tminIdx)));

    % Setup colors.
    startColor = [79., 0., 128.]/255.;
    endColor = [255., 143., 255.]/255.;
    r = linspace(startColor(1), endColor(1), length(times(1:tminIdx)));
    g = linspace(startColor(2), endColor(2), length(times(1:tminIdx)));
    b = linspace(startColor(3), endColor(3), length(times(1:tminIdx)));
    xcurr = human.x;

    % Record state.
    X(1) = xcurr(1);
    Y(1) = xcurr(2);
    PG(1) = xcurr(3);

    % Plot first point.
    color = [r(1), g(1), b(1)];
    plot3(xcurr(1), xcurr(2), xcurr(3), '-o', 'color', color, 'markeredgecolor', color, 'markerfacecolor', color);
    for t=2:length(times(1:tminIdx))
        % Apply control to dynamics.
        human.updateState(uopt{t}, times(t), human.x);
        xcurr = human.x;

        % record state.
        X(t) = xcurr(1);
        Y(t) = xcurr(2);
        PG(t) = xcurr(3);

        % Plot point and connection between pts.
        color = [r(t), g(t), b(t)];
        p = plot3([X(t-1), X(t)], [Y(t-1), Y(t)], [PG(t-1), PG(t)], '-o', ...
                'Color', color, ...
                'markeredgecolor', color, ...
                'markerfacecolor', color);
        p.LineWidth = 2;
        
        % add timestamps
        txt = num2str(times(t));
        t = text(X(t)+0.05, Y(t)+0.05, PG(t)+0.05, txt);
        t.Color = color;
    end

    % Plot goals (red is ground truth, grey is other goal).
    plot3(goals{trueGoalIdx}(1), goals{trueGoalIdx}(2), 0.5, '-o', ...
                'Color', 'r', ...
                'markeredgecolor', 'r', ...
                'markerfacecolor', 'r');
    g1Txt = strcat('g', num2str(trueGoalIdx));
    t1 = text(goals{trueGoalIdx}(1), goals{trueGoalIdx}(2), 0.55, g1Txt);
    t1.FontSize = 12;
    t1.Color = 'r';

    if trueGoalIdx == 2
        otherGoalIdx = 1; 
    else
        otherGoalIdx = 2;
    end
    plot3(goals{otherGoalIdx}(1), goals{otherGoalIdx}(2), 0.5, '-o', ...
                'Color', 'k', ...
                'markeredgecolor', 'k', ...
                'markerfacecolor', 'k');
            
    g2Txt = strcat('g', num2str(otherGoalIdx));
    t2 = text(goals{otherGoalIdx}(1), goals{otherGoalIdx}(2), 0.55, g2Txt);
    t2.FontSize = 12;
    t2.Color = 'k';
    
    grid on;
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    zlim([grid_min(3), grid_max(3)]);
    
    xlabel('x');
    ylabel('y');
    zlabel('P(g = g1)');
end