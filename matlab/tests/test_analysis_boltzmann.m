clc
clf
clear all

%% Grid
grid_min = [-5; -5; -0.1];  % Lower corner of computation domain
grid_max = [5; 5; 1.1];     % Upper corner of computation domain
N = [11; 11; 11];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Create human dynamical system
% Case 1: Human's beta* = 1 (Boltzmann optimal)
%         Robot's Prior(beta = 0) = 0.1 (robot thinks human is random)
%          
%         Analysis 1.1: min_u computes the min time it takes robot to
%                       realize human is actually optimal
%         Analysis 1.2: max_u computes the max time it takes robot to
%                       realize human is actually optimal
%
% Case 2: Human's beta* = 0 (random)
%         Robot's Prior(beta = 1) = 0.9 (robot thinks human is optimal)
%
%         Analysis 2.1: min_u computes the min time it takes robot to
%                       realize human is actually random
%         Analysis 2.2: max_u computes the max time it takes robot to
%                       realize human is actually random

% Velocity
v = 0.2;

% Control bounds
uRange = [-pi+1e-2; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 1;

% Number of discrete controls
numCtrls = 50;

% Threshold to determine likely controls
uThresh = 0;

% Are we using dynamic of static beta model?
betaModel = 'static';

% Human's goal location.
theta = [2, 2];

% Timestep in discretized dynamics (for Q-function computation).
delta_t = 0.1;

% We have no dynamic beta parameters
extraArgs = [];

% Tolerance for how "sufficiently high" the probability needs to be
tol = 0.2;

% ---- Setup for Case 1 ---- %
% start with high prior on beta=0, but true is beta=1 (boltzmann)
Pbeta1 = 0.8; 
trueBeta = 1;
centerPBeta = 1;

% For Analysis 1.1 (min time)
uMode = 'min';

% For Analysis 1.2 (max time)
%uMode = 'max';
% -------------------------- %

% % ---- Setup for Case 2 ---- %
% % start with high prior on beta=1, but true is beta = 0 (rand)
% Pbeta1 = 0.9; 
% trueBeta = 0;
% centerBeta = 0;
% 
% % For Analysis 2.1 (min time)
% uMode = 'min';
% 
% % For Analysis 2.2 (max time)
% % Note: set should not grow because the random human can always choose the 
% %       optimal action (u*=0) that is most likely under the Gaussian human,
% %       thereby fooling you into thinking that he's Gaussian. 
% %uMode = 'max';  
% % -------------------------- %

% Setup dynamical system
x0 = [-3; 0; Pbeta1];
human = Boltzmann1or0Human(x0, v, trueBeta, uRange, gamma, theta, ...
    delta_t, uThresh, numCtrls, betaModel, extraArgs);

%% Setup target set
% Target set is centered at the true beta value
xyoffset = 0.1;
poffset = 0.01;
center = [0; 0; centerPBeta];
widths = [(grid_max(1) - grid_min(1)) - xyoffset; ...
          (grid_max(2) - grid_min(2)) - xyoffset; 
          tol - poffset];

% center = [theta(1); theta(2); centerPBeta];
% widths = [1; ...
%           1; 
%           tol - poffset];

% center = [theta(1); theta(2); 0.5];
% widths = [1; ...
%           1; 
%           (grid_max(1) - grid_min(1)) + tol - poffset];
data0 = shapeRectangleByCenter(g, center, widths);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% PLOTS FOR DEBUGGING.
% human.plotPUGivenXBeta(g);

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
schemeData.hamFunc = @boltzmann1or0Human_ham;
schemeData.partialFunc = @boltzmann1or0Human_partial;

%% Setup value function computation params.
% HJIextraArgs.visualize = true; %show plot
HJIextraArgs.quiet = false;
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [grid_min(1) grid_max(1) ...
                                   grid_min(2) grid_max(2) ...
                                   -0.1 1.1];
HJIextraArgs.visualize.xTitle = '$p^x$';
HJIextraArgs.visualize.yTitle = '$p^y$';
HJIextraArgs.visualize.zTitle = '$P(\beta = 1) [Boltz]$';
HJIextraArgs.visualize.fontSize = 15;
HJIextraArgs.stopInit = x0;
%HJIextraArgs.visualize.camlightPosition = [0 0 0];

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid
HJIextraArgs.ignoreBoundary = 0; 

% minWith = 'set';
minWith = 'zero';
%minWith = 'minVwithL';

%% Debugging
fprintf('------ Prediction Analysis Setup -------\n');
fprintf('   true human beta: %d\n', trueBeta);
fprintf('   [x, y, P(beta=1)]: [%d, %d, %d]\n', x0(1), x0(2), x0(3));
fprintf('   uMode: %s\n', uMode);
fprintf('--------------------------------------\n');

%% Solve it!
[valueFuns, times, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);
tminIdx = length(times);

%% Grab the min/max time. 
tmin = -Inf;
for t=1:length(times)
    v = eval_u(g, valueFuns(:,:,:,t), x0); 
    if v <= 0
        tmin = (t-1)*dt;
        break;
    end
end

fprintf("Minimum time it takes to realize beta=%d is %f\n", trueBeta, tmin);

%% Get the sequence of optimal controls.
uopt = GetOptControls(x0, g, valueFuns, times(1:tminIdx), human, uMode);

%% Simulate and visualize optimal control and states.
goalSetRad = 1;
plotStateTraj(uopt, human, times, tminIdx, theta, ...
    grid_min, grid_max, goalSetRad, dt);

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

%% Plots the states that result from optimal control sequence.
function plotStateTraj(uopt, human, times, tminIdx, goal, ...
    grid_min, grid_max, goalSetRad, dt)
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
        human.updateState(uopt{t}, dt, human.x);
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
%         txt = num2str(times(t));
%         t = text(X(t)+0.05, Y(t)+0.05, PG(t)+0.05, txt);
%         t.Color = color;
    end


    %add timestamps
    txt = strcat('t=', num2str(times(t)), ', p=', num2str(PG(t)));
    t = text(X(t)+0.2, Y(t)+0.05, PG(t)+0.05, txt);
    t.Color = color;

    % Plot goal
    rectangle('Position',[goal(1)-goalSetRad ...
                          goal(2)-goalSetRad ...
                          goalSetRad*2 ...
                          goalSetRad*2],...
                          'Curvature',1, ...
                          'FaceColor',[1, 0.67, 0.67],...
                          'EdgeColor',[1, 0.67, 0.67],...
                          'LineWidth',1);
    plot3(goal(1), goal(2), 0.5, '-o', ...
                'Color', 'r', ...
                'markeredgecolor', 'r', ...
                'markerfacecolor', 'r');
    g1Txt = 'g';
    t1 = text(goal(1), goal(2), 0.55, g1Txt);
    t1.FontSize = 12;
    t1.Color = 'r';
    
    grid on;
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    zlim([grid_min(3), grid_max(3)]);
    
    xlabel('x');
    ylabel('y');
    zlabel('P(g = g1)');
end