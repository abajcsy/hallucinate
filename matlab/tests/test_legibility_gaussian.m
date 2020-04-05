clc
clf
close all
clear all

%% Grid
grid_min = [-10; -10; -0.1];  % Lower corner of computation domain
grid_max = [10; 10; 1.1];     % Upper corner of computation domain
N = [31; 31; 31];           % Number of grid points per dimension
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
piTol = 1e-2;
uRange = [-pi+piTol; pi-piTol]; 

% gamma in continuous-time P(beta = 0) dynamics
gamma = 0.1;

% Number of discrete controls
numCtrls = 8;

% Variance in normal distributions
sigma = pi/8;

% Known human goal locations. 
%goals = {[5,-5], [5,5]};
goals = {[2,-2], [2,2]};
%goals = {[2,-0.5], [2,0.5]};
%goals = {[2, 0], [4, 0]};

% Threshold to determine likely controls
uThresh = 0.01;  

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
trueGoalIdx = 1; %2;
goalSetRad = 0.5;
centerPgoal1 = 1; %0;

% For Analysis 1.1 & 2.1 (min time)
uMode = 'min';

% For Analysis 1.2 & 2.2 (max time)
%uMode = 'max';

% ---- Plotting info --- %
extraPltArgs.compType = 'conf';
%extraPltArgs.compType = 'goal';
%extraPltArgs.compType = 'conf_and_goal';
extraPltArgs.uThresh = uThresh;
extraPltArgs.uMode = uMode;
% ---- Plotting info --- %

% Setup dynamical system
x0 = [0; 0; Pgoal1];
human = GaussianG1orG2Human(x0, v, trueGoalIdx, goalSetRad, uRange, gamma, goals, ...
                        sigma, uThresh, numCtrls, betaModel, piTol, extraArgs);

%% Setup target set
% Target set is centered at the true beta value
xyoffset = 0.1;
poffset = 0.01;

if strcmp(extraPltArgs.compType, 'conf')
    % Target = anywhere in X/Y and high confidence. 
    center = [0; 0; centerPgoal1];
    widths = [(grid_max(1) - grid_min(1)) - xyoffset; ...
              (grid_max(2) - grid_min(2)) - xyoffset; 
              tol - poffset];
    data0 = shapeRectangleByCenter(g, center, widths);
elseif strcmp(extraPltArgs.compType, 'conf_and_goal')
    %Target = goal set and high confidence.  
    center = [goals{trueGoalIdx}(1); goals{trueGoalIdx}(2); centerPgoal1];
    widths = [goalSetRad*2; ...
              goalSetRad*2; 
              0.2];
    data0 = shapeRectangleByCenter(g, center, widths);
elseif strcmp(extraPltArgs.compType, 'goal')
    %Target = goal set and any confidence.
    data0 = shapeCylinder(g, 3, ...
        [goals{trueGoalIdx}(1); goals{trueGoalIdx}(2); centerPgoal1], goalSetRad);
else
    error("No supported computation type!");
end

%% Let the human have access to the grid for debugging.
human.setGrid(g);

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
schemeData.accuracy = 'high'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'backward';
schemeData.hamFunc = @gaussianG1orG2Human_ham;
schemeData.partialFunc = @gaussianG1orG2Human_partial;

%% Setup value function computation params.
HJIextraArgs.quiet = false;
HJIextraArgs.visualize.valueSet = 1;
HJIextraArgs.visualize.initialValueSet = 0;
HJIextraArgs.visualize.figNum = 1; %set figure number
HJIextraArgs.visualize.deleteLastPlot = true; %delete previous plot as you update
HJIextraArgs.visualize.viewGrid = true;
HJIextraArgs.visualize.viewAxis = [grid_min(1) grid_max(1) ...
                                   grid_min(2) grid_max(2) ...
                                   -0.1 1.1];
HJIextraArgs.visualize.xTitle = '$x$';
HJIextraArgs.visualize.yTitle = '$y$';
HJIextraArgs.visualize.zTitle = '$P(goal = 1)$';
HJIextraArgs.visualize.fontSize = 15;
HJIextraArgs.stopInit = x0;
HJIextraArgs.visualize.plotColorVS = [0.9,0.9,0.9];

% % Uncomment to see a 2D slice
% HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
% HJIextraArgs.visualize.plotData.projpt = [0.5]; 
% HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

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
% Note: The LevelSetToolbox can only solve INITIAL-value PDEs. However, 
% we are solving for a BRS, which is a FINAL-vablue PDE:
%          
%       dV/dt + min_u \grad V(x,t) * f(x,u) = 0
%       V(x,T) = l(x)
%
% Because we are computing a BRS, the order of the value functions is:
%
%       valueFuns(1)    = V(x,T)
%       valueFuns(2)    = V(x, T-dt)
%       ...
%       valueFuns(end)  = V(x, 0)
%
% The values in the times variable go from:
%
%       times(1)    = 0
%       times(2)    = dt
%       times(3)    = 2*dt
%       ...
%       times(end)  = T
%
% However, the "real" time that corresponds to the V's should have the time
% flipped. 
[valueFuns, times, ~] = ...
  HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);

% Flip the order of the value functions to match the time. This way we
% start from the beginning of time: V(x,0) and go to V(x,T), just like
% time.
valueFuns = flip(valueFuns,4);

% %% ------- DEBUGGING CODE. --------- %%
% figure(3);
% clf;
% timeIdx = 10;
% visSetIm(g, valueFuns(:,:,:,timeIdx));
% title(strcat("Value Function at t=", num2str(times(timeIdx))));
% xlim([grid_min(1) grid_max(1)]);
% ylim([grid_min(2) grid_max(2)]);
% zlim([0,1]);
% grid on;
% %------- DEBUGGING CODE. --------- %

%% Grab the min time. 
tmin = times(end);
tminIdx = length(times);

%% Get the optimal trajectory.
[traj, traj_tau] = computeOptTraj(g, valueFuns, times, human, extraArgs);

fprintf("Minimum time it takes to realize goal=%d is %f\n", ...
    trueGoalIdx, traj_tau(end));
        
%% Plot the trajectory.
plotTraj(traj, traj_tau, goals, trueGoalIdx, ...
    grid_min, grid_max, goalSetRad, extraPltArgs);

% save('uopt_09thresh.mat', 'uopt', 'human', 'times', 'tminIdx', 'goals', 'trueGoalIdx', 'grid_min', 'grid_max', 'goalSetRad', 'dt');

%% Plots the state trajectory.
function plotTraj(traj, traj_tau, goals, trueGoalIdx, ...
    grid_min, grid_max, goalSetRad, extraPltArgs)
    figure(2);
    hold on
    
    % Goal colors.
    g1Color = 'r';
    g2Color = [38., 138., 240.]/255.; %'b';

    % Setup colors.
    %startColor = [79., 0., 128.]/255.;
    %endColor = [255., 143., 255.]/255.;
    startColor = [97., 76., 76.]/255.;
    endColorRed = [255., 0., 0.]/255.;
    endColorBlue = [38., 138., 240.]/255.;
    
    if trueGoalIdx == 1
        endColor = endColorRed;
        finaloffset = -0.3;
        initoffset = 0.3;
    else
        endColor = endColorBlue;
        finaloffset = 0.3;
        initoffset = -0.3;
    end
    
    r = linspace(startColor(1), endColor(1), length(traj_tau));
    g = linspace(startColor(2), endColor(2), length(traj_tau));
    b = linspace(startColor(3), endColor(3), length(traj_tau));

    % Record state.
    % Plot first point.
    color = [r(1), g(1), b(1)];
    xcurr = traj(1:3, 1);
    plot3(xcurr(1), xcurr(2), xcurr(3), '-o', 'color', color, ...
        'markeredgecolor', color, 'markerfacecolor', color);
    
    % Add first timestamp.
    txt = strcat('t=', num2str(traj_tau(1)), ', p=', num2str(xcurr(3)));
    tp = text(xcurr(1)+0.05, xcurr(2)+initoffset, xcurr(3)+0.05, txt);
    tp.Color = color;
    for t=2:length(traj_tau)
        xprev = traj(1:3, t-1);
        xcurr = traj(1:3, t);
        % Plot point and connection between pts.
        color = [r(t), g(t), b(t)];
        p = plot3([xprev(1), xcurr(1)], [xprev(2), xcurr(2)], [xprev(3), xcurr(3)], '-o', ...
                'Color', color, ...
                'markeredgecolor', color, ...
                'markerfacecolor', color);
        p.LineWidth = 2;
    end
    
    xcurr = traj(1:3, end);
    % add final info.
    txt = strcat('t=', num2str(traj_tau(end)), ', p=', num2str(xcurr(3)));
    tp = text(xcurr(1), xcurr(2)+finaloffset, xcurr(3)+0.1, txt);
    tp.Color = color;

    % Plot GOAL 1.
    if strcmp(extraPltArgs.compType, 'conf_and_goal') || ...
            strcmp(extraPltArgs.compType, 'goal')
        
        rectangle('Position',[goals{1}(1)-goalSetRad ...
                              goals{1}(2)-goalSetRad ...
                              goalSetRad*2 ...
                              goalSetRad*2],...
                              'Curvature',1, ...
                              'FaceColor',[1, 0.67, 0.67],...
                              'EdgeColor',[1, 0.67, 0.67],...
                              'LineWidth',1);
    end
    plot3(goals{1}(1), goals{1}(2), 0.5, '-o', ...
                'Color', g1Color, ...
                'markeredgecolor', g1Color, ...
                'markerfacecolor', g1Color);
    g1Txt = 'g1';
    t1 = text(goals{1}(1)+0.3, goals{1}(2), 0.55, g1Txt);
    t1.FontSize = 12;
    t1.Color = g1Color;

    % Plot GOAL 2
    if strcmp(extraPltArgs.compType, 'conf_and_goal') || ...
            strcmp(extraPltArgs.compType, 'goal')
        rectangle('Position',[goals{2}(1)-goalSetRad ...
                  goals{2}(2)-goalSetRad ...
                  goalSetRad*2 ...
                  goalSetRad*2],...
                  'Curvature',1, ...
                  'FaceColor',[0.7098, 0.8980, 1.0000],...
                  'EdgeColor',[0.7098, 0.8980, 1.0000],...
                  'LineWidth',1);
    end
    plot3(goals{2}(1), goals{2}(2), 0.5, '-o', ...
                'Color', g2Color, ...
                'markeredgecolor', g2Color, ...
                'markerfacecolor', g2Color);
            
    g2Txt = 'g2';
    t2 = text(goals{2}(1)+0.3, goals{2}(2), 0.55, g2Txt);
    t2.FontSize = 12;
    t2.Color = g2Color;
    
    grid on;
    xticks([grid_min(1), -3, -2, -1, 0, 1, 2, 3, grid_max(1)]);
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    zlim([grid_min(3), grid_max(3)]);
    
    xlabel('x');
    ylabel('y');
    zlabel('P(g = g1)');
    set(gcf,'Position',[100 100 500 500]);
    
    % Saving functionality.
    repo = what('hallucinate');
    filename = strcat('g', num2str(trueGoalIdx), '_uthr' , ...
        num2str(extraPltArgs.uThresh), '_', extraPltArgs.compType, ...
        '_', extraPltArgs.uMode, '.png');
    saveas(gcf, strcat(repo.path, '/legibility_imgs/', filename));
    
    view(17, 14);
    filename = strcat('g', num2str(trueGoalIdx), '_uthr' , ...
        num2str(extraPltArgs.uThresh), '_', extraPltArgs.compType, ...
        '_', extraPltArgs.uMode, '_view2.png');
    saveas(gcf, strcat(repo.path, '/legibility_imgs/', filename));
end

