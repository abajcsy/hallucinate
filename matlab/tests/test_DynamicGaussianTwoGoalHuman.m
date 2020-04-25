clc
clear all

%% Create human dynamical system

% Velocity
v = 0.6;

% Control bounds
uRange = [-pi; pi];

% gamma in continuous-time P(beta = 0) dynamics
gamma = 0.01;

% Number of discrete controls
numCtrls = 31;

% Threshold to determine likely controls
uThresh = 0.01;

% Variance in normal distributions
sigma = pi/4;

% Known human goal locations. 
goals = {[2, -2], [2,2]}; 

% Are we using dynamic of static beta model?
%betaModel = 'static';
betaModel = 'dynamic';

% Dynamic beta parameters
%extraArgs = [];
extraArgs.DeltaB0 = 0.5; % beta = 0 <==> goal = 1
extraArgs.alpha = 0.8;

% Setup dynamical system
Pgoal1 = 0.9; 
x0 = [0; 0; Pgoal1];
human = GaussianTwoGoalHuman(x0, v, uRange, gamma, goals, sigma, uThresh, numCtrls, ...
    betaModel, extraArgs);

human.debugMode = false;
human.percentileThresh = 0.99;
% human.percentileThresh = 0.10;
human.usePercentileThresh = true;

%% Grid
grid_min = [-4; -4; -0.1];  % Lower corner of computation domain
grid_max = [4; 4; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Let the human have access to the grid for debugging.
human.setGrid(g);

%% Pre-compute the optimal control over the entire state-space.
human.computeUOptGoals(g.xs);

%% Pre-compute the likely controls and dynamics over the entire state-space.
human.computeUAndXDot(g.xs);

%% target set
% Findings: (1) Increasing R reduces the gap between the two extremes so it 
% could just be numerical issues. (2) Increasing the grid resolution for
% R=0.1 also reduces the gap between the two extremes so that is other
% evidence that it is just numerical error.
R = 0.1;
data0 = shapeSphere(g, x0, R);

%% time vector
t0 = 0;
tMax = 1.8;
dt = 0.1667; %0.05;
tau = t0:dt:tMax;
uMode = 'max';

%% Pack problem parameters
% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = human;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'forward';
schemeData.hamFunc = @gaussianTwoGoalHuman_ham;
schemeData.partialFunc = @gaussianTwoGoalHuman_partial;

%% Compute value function
% HJIextraArgs.visualize = true; %show plot
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = false; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [grid_min(1) grid_max(1) ...
                                   grid_min(2) grid_max(2) ...
                                   -0.1 1.1];
HJIextraArgs.visualize.xTitle = '$p^x$';
HJIextraArgs.visualize.yTitle = '$p^y$';
HJIextraArgs.visualize.zTitle = '$P(goal_1)$';
HJIextraArgs.visualize.fontSize = 15;
%HJIextraArgs.visualize.camlightPosition = [0 0 0];

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid
HJIextraArgs.ignoreBoundary = 0; 

%uncomment if you want to see a 2D slice
HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
HJIextraArgs.visualize.plotData.projpt = {'min'}; %project pt
HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

%HJIextraArgs.targets = data0;

minWith = 'set';
%minWith = 'zero';
%minWith = 'minVwithL';

tStart = tic;
[data, tau2, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);
tEnd = toc(tStart);
fprintf(strcat("[Reachability] Dynamic param pred time: ", num2str(tEnd), " s\n"));

%% Plotting!
figure(3);
hold on

% Color setup.
start_color = [97, 17, 44]/255.;
end_color = [237, 19, 92]/255.;
red = linspace(start_color(1), end_color(1), length(tau2));
green = linspace(start_color(2), end_color(2), length(tau2));
blue = linspace(start_color(3), end_color(3), length(tau2));

dimsToRemove = [0 0 1];
for t=1:length(tau2)
    
    titleString = strcat('Reachability: Dynamic Param. HMM alpha=', ...
        num2str(extraArgs.alpha), ', t=', num2str(tau2(t)), 's');
    title(titleString);
    
    [g2d, data2D] = proj(g, data(:,:,:,t), dimsToRemove, 'min');
    
    figure(3);
    hold on;
    
    %Plot prediction contour.
    [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
    c.LineWidth = 2;
    c.EdgeColor = [red(t), green(t), blue(t)];
    
    % Plot goals.
    figure(3);
    scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
    scatter(goals{2}(1), goals{2}(2), 100, 'b', 'filled');
    
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    set(gcf,'Position',[100 100 700 700]);
    set(gcf,'color','w');
    whitebg('k');
    grid on
    %pause(0.1);
end

