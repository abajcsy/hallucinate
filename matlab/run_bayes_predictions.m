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

%% (Bayesian) Predictor params.

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

% (reachablity) Threshold to determine likely controls -- used for
% plotting.
uThresh = 0.1; 

% Variance on Gaussian observation model.
sigma1 = pi/4;
sigma2 = pi/4;

% Set the prior over goal 1 and goal 2.
prior = [0.5, 0.5];
Pgoal1 = prior(1); 

% Grid cell size.
r = 0.1;

type_of_pred = 'static';
%type_of_pred = 'dynamic';

% Create the predictor.
if strcmp(type_of_pred, 'static')
    predictor = LatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
        grid_min, grid_max, r);
elseif strcmp(type_of_pred, 'dynamic') 
    % Discrete distribution over two goals: \delta{1,2}.
    Delta = [0.5, 0.5];

    % Probability of staying in same goal (in HMM model).
    hmmParam = 0.8;

    % Create the predictor. 
    predictor = DynamicLatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
        grid_min, grid_max, r, hmmParam, Delta);
else
    error("Illegal type of pred.");
end

% Initial human state (in m).
x0 = [0,0];
v = 0.6;
dt = r / v;

% Prediction horizon. 
% gString = createGrid(gridMin, gridMax, gridDims);
% dt = gString.dx(1)/v;
T = 1.8;        %3              % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)
discrete_times = 0:1:H+1;       % (used for data storing and plotting)

%% Simulation params. 

% Number of simulation steps. (real sim time = T*dt)
simT = 36;

% State. 
xcurr = x0;

% Color setup.
start_color = [0, 77, 73]/255.;
end_color = [0, 232, 220]/255.;
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
    preds = predictor.predict(xcurr, H);
    all_preds{end+1} = preds; 
    all_pred_times{end+1} = discrete_times*dt;

    %% Plotting. 
    if plotData
        % Plot the predictions.
        figure(2);
        hold on;
        for pt=1:length(preds)
            curr_pred = preds{pt};
            [opt_eps, P, X, Y] = ...
                compute_likely_states(curr_pred, predictor, uThresh);

            % Plot prediction contour.
            [M, c] = contour(X, Y, P, [1, 1]);
            c.LineWidth = 2;
            c.EdgeColor = [red(t), green(t), blue(t)];
        end

        % Plot state of human.
        scatter(xcurr(1), xcurr(2), 60, [hcolor(t), hcolor(t), hcolor(t)], 'filled');

        % Plot goals.
        figure(2);
        scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
        scatter(goals{2}(1), goals{2}(2), 100, 'b', 'filled');

        scatter(trueGoal{1}(1), trueGoal{1}(2), 100, 'y', 'filled');

        xlim([grid_min(1), grid_max(1)]);
        ylim([grid_min(2), grid_max(2)]);
        set(gcf,'Position',[100 100 700 700]);
        set(gcf,'color','w');
        whitebg('k');
        grid on
        %pause(dt);

        % Plot posterior. 
        figure(3);
        hold on
        plot(times, posteriors, 'r-o', 'LineWidth', 3);
        xlabel("$time$", 'Interpreter', 'Latex');
        ylabel("$P(g_1)$", 'Interpreter', 'Latex');
        ylim([0,1]);
        grid on
    end
    
    %% Update. 
    % Forward simulate human. 
    [xnext, ucurr] = human.simulate(xcurr, t*dt, dt);
    xcurr = xnext;
    
    % Update Posterior.
    [i,j] = predictor.realToSim(xcurr);
    s = [i,j];
    pu_goal1 = predictor.Pu_given_x_g_normalized(ucurr, s, 1);
    pu_goal2 = predictor.Pu_given_x_g_normalized(ucurr, s, 2);
    Pgoal1 = (pu_goal1*Pgoal1)/(pu_goal1*Pgoal1 + pu_goal2*(1-Pgoal1));
    
    % Update predictor info.
    predictor.prior = [Pgoal1, (1-Pgoal1)];
end

if saveData
    repo = what('hallucinate');
    filename = strcat('bayesian_', type_of_pred, '_pg1', num2str(prior(1)), '.mat');
    save(strcat(repo.path, '/ral_data/', filename), ...
        'human_states', 'times', 'posteriors', 'all_preds', ...
        'all_pred_times', 'predictor', 'prior');
end

%% Grab all the likely-enough predicted states.
function [opt_eps, P, X, Y] = compute_likely_states(preds, predictor, ...
    delta_reachability)
    
    % Grid for Bayesian prediction
    [X, Y] = predictor.getLatticeMeshgrid();
    
    valid_indices = find(preds > 0);
    valid_data = preds(valid_indices);
    sorted_valid_data = sort(valid_data, 'descend');
    eps_index = find(cumsum(sorted_valid_data) > (1 - delta_reachability), 1, 'first');
    
    if isempty(eps_index)
        % if we can't find likely enough states, then we should take max. 
        eps_index = find(max(cumsum(sorted_valid_data)), 1, 'first');
    end
    
    opt_eps = sorted_valid_data(eps_index);
    
    % Compute the optimal predictions
    P = zeros(size(X));
    for r = 1:predictor.rows
        for c = 1:predictor.cols
            linIdx = sub2ind([predictor.rows,predictor.cols],r,c);
            if delta_reachability > 0.0
                P(r, c) = 1*(preds(linIdx) >= opt_eps) + 0*(preds(linIdx) < opt_eps);
            else
                P(r, c) = 1*(preds(linIdx) > 0.0);
            end
        end
    end
end


