clf
%% Random seed.
rng(3);

%% Simulated human params.

% Velocity
v = 0.6;

% Control bounds
uRange = [-pi; pi];

% Standard dev in normal distributions
%human_sigma = pi/4;
human_sigma = 0.0;

% Known human goal locations. 
goals = {[2, 2], [2, -2]}; 

% True goal that human is going to.
trueIdx = 1;
trueGoal = goals(trueIdx);
%trueGoal = {[3, 0]};

% Initial state of human.
x0 = [0; 0];

% Create simulated human.
human = GaussianGoalHuman(x0, v, human_sigma, uRange, trueGoal);

%% (Reachability) Predictor params.

% Grid
grid_min = [-4; -4; -0.1];  % Lower corner of computation domain
grid_max = [4; 4; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

% gamma in continuous-time P(beta = 0) dynamics
gamma = 0.5;

% Number of discrete controls
numCtrls = 31;

% Threshold to determine likely controls
uThresh = 0.005; 

% Are we using dynamic of static beta model?
type_of_pred = 'static';
%type_of_pred = 'dynamic';

if strcmp(type_of_pred, 'static')
    betaModel = 'static';
    % (No) Dynamic beta parameters
    extraArgs = [];
elseif strcmp(type_of_pred, 'dynamic')
    betaModel = 'dynamic';
    % Dynamic beta parameters
    extraArgs.DeltaB0 = 0.5; % beta = 0 <==> goal = 1
    extraArgs.alpha = 0.8;
end

% Prediction sigma.
sigma = pi/4;

% Prior. 
prior = [0.5,0.5];
Pgoal1 = prior(1); 

% Setup dynamical system
z0 = [x0; Pgoal1];
predictor = GaussianTwoGoalHuman(z0, v, uRange, gamma, goals, ...
    sigma, uThresh, numCtrls, ...
    betaModel, extraArgs);
predictor.debugMode = false;
predictor.usePercentileThresh = false;

% Let the human have access to the grid for debugging.
predictor.setGrid(g);

% Pre-compute the optimal control over the entire state-space.
fprintf("Pre-computing opt control for each goal over all states...\n");
predictor.computeUOptGoals(g.xs);

% Pre-compute the likely controls and dynamics over the entire state-space.
fprintf("Pre-computing dynamics for each state and control...\n");
predictor.computeUAndXDot(g.xs);

% Target set.
R = 0.1;
data0 = shapeSphere(g, z0, R);

% Time vector
t0 = 0;
tMax = 3; %1.8;
dt = 0.1667; %0.05;
tau = t0:dt:tMax;
uMode = 'max';

% Put grid and dynamic systems into schemeData
schemeData.grid = g;
schemeData.dynSys = predictor;
schemeData.accuracy = 'medium'; %set accuracy
schemeData.uMode = uMode;
schemeData.tMode = 'forward';
schemeData.hamFunc = @gaussianTwoGoalHuman_ham;
schemeData.partialFunc = @gaussianTwoGoalHuman_partial;

% Setup extraArgs.
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

% % uncomment if you want to see a 2D slice
% HJIextraArgs.visualize.plotData.plotDims = [1 1 0]; %plot x, y
% HJIextraArgs.visualize.plotData.projpt = {'min'}; %project pt
% HJIextraArgs.visualize.viewAngle = [0,90]; % view 2D

% since we have a finite compute grid, we may not want to 
% trust values near the boundary of grid.
HJIextraArgs.ignoreBoundary = 0; 

% Compute FRS.
minWith = 'set';

%% Simulation params. 

% Number of simulation steps. (real sim time = T*dt)
simT = 0; %36;

dimsToRemove = [0 0 1];
xcurr = x0;

% Color setup.
start_color = [97, 17, 90]/255.;
end_color = [235, 19, 216]/255.;
red = linspace(start_color(1), end_color(1), simT+1);
green = linspace(start_color(2), end_color(2), simT+1);
blue = linspace(start_color(3), end_color(3), simT+1);
hcolor = linspace(0.3, 1, simT+1);

% Saving and plotting variables.
human_states = {};
times = [];
posteriors = [];
all_preds = {};
all_pred_times = {};

%% Save data? Plot?
saveData = true;
plotData = true;

%% Simulate.
hold on;
for t=1:simT+1
    % Store current state.
    human_states{end+1} = xcurr;
    
    % Store time and posterior.
    times(end+1) = (t-1)*dt;
    posteriors(end+1) = Pgoal1;
    
    % Predict the human.
    figure(1);
    [preds, pred_times, ~] = ...
        HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);
    all_preds{end+1} = preds; 
    all_pred_times{end+1} = pred_times; 
    
    %% Plotting. 
    if plotData
        % Plot the predictions.
        figure(2);
        hold on;
        for pt=1:length(pred_times)
            [g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
            [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
            c.LineWidth = 2;
            c.EdgeColor = [red(t), green(t), blue(t)];
        end

        % Plot state of human.
        scatter(xcurr(1), xcurr(2), 60, [hcolor(t), hcolor(t), hcolor(t)], 'filled');

        % Plot goals.
        figure(2);
        scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
        scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');

        scatter(trueGoal{1}(1), trueGoal{1}(2), 100, 'y', 'filled');

        xlim([grid_min(1), grid_max(1)]);
        ylim([grid_min(2), grid_max(2)]);
        set(gcf,'Position',[100 100 700 700]);
        set(gcf,'color','w');
        whitebg('w');
        grid on

        % Plot posterior. 
        figure(3);
        hold on
        plot(times, posteriors, 'r-o', 'LineWidth', 3);
        xlabel("$time$", 'Interpreter', 'Latex');
        ylabel("$P(g_1)$", 'Interpreter', 'Latex');
        ylim([0,1]);
        grid on
    end
    
    %% Forward simulate human. 
    [xnext, ucurr] = human.simulate(xcurr, t*dt, dt);
    xcurr = xnext;
    
    %% Update Posterior.
    xcurr_cell = {xcurr(1), xcurr(2)}; 
    pu_goal1 = predictor.PuGivenGoal_normalized(ucurr, xcurr_cell, 1);
    pu_goal2 = predictor.PuGivenGoal_normalized(ucurr, xcurr_cell, 2);
    Pgoal1 = (pu_goal1*Pgoal1)/(pu_goal1*Pgoal1 + pu_goal2*(1-Pgoal1));
    
    %% Update initial set. 
    z0 = [xcurr(1); xcurr(2); Pgoal1];
    data0 = shapeSphere(g, z0, R);
    
    %% Update predictor info.
    predictor.x = z0;
    schemeData.dynSys = predictor;
end

if saveData
    repo = what('hallucinate');
    filename = strcat('reach_', type_of_pred, '_pg1', num2str(prior(1)), ...
        '_gamma', num2str(gamma), '_delta', num2str(uThresh), '_1pred.mat');
    save(strcat(repo.path, '/ral_data/', filename), ...
        'human_states', 'times', 'posteriors', 'all_preds', ...
        'all_pred_times', 'predictor', 'prior', 'uThresh', 'g');
end

