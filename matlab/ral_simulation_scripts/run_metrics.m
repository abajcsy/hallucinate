clear all
close all

%% Load data.
repo = what('hallucinate');
filename = 'sigma_init_cond.mat';
load(strcat(repo.path, '/ral_data_revise/human_traj_data/', filename));

%% (Bayesian) Predictor params.

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

% (reachablity) Threshold to determine likely controls -- used for
% plotting.
uThresh = 0.03; 

% Variance on Gaussian observation model.
sigma1 = pi/4;
sigma2 = pi/4;

% Set the prior over goal 1 and goal 2.
prior = [0.5, 0.5];

% Grid cell size.
r = 0.1;

% Velocity
v = 0.6;

% Known human goal locations. 
goals = {[2; 2], [2; -2]}; 

% Create the predictor.
predictor = LatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
    grid_min, grid_max, r);

dt = r / v;
T = 2;                          % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)
discrete_times = 0:1:H+1;       % (used for data storing and plotting)

all_con = [];
all_acc = [];
%% Run metrics!
for sigma = human_sigmas
    traj_map = all_trajs(num2str(sigma));
    for init_state = human_init_conds
        x0 = init_state{1};
        human_traj = traj_map(num2str(x0));
        
        % Main simulation loop!
        for t=1:simT
            xcurr = human_traj{t};
            
            % Predict the human.
            preds = predictor.predict(xcurr, H);
            
            % Compute % overconservative. 
            p_con = conservativism(preds, human_traj, v, dt, predictor, uThresh);
            
            % Compute % accurate. 
            p_acc = accuracy(preds, xcurr, t, human_traj, predictor, uThresh);
            
            fprintf('for sigma=%f, init=(%f,%f): \n', sigma, x0(1), x0(2));
            fprintf('    CON = %f\n', p_con);
            fprintf('    ACC = %f\n', p_acc);
            all_con(end+1) = p_con;
            all_acc(end+1) = p_acc;
        end
    end
end

% Return the average over all accuracy and conservativism!
avg_con = mean(all_con);
avg_acc = mean(all_acc);

fprintf('AVG CON = %f\n', avg_con);
fprintf('AVG ACC = %f\n', avg_acc);

%% Returns percent of the full FRS that the predictions occupy.
% The lower the better.
function p_con = conservativism(preds, human_traj, v, dt, predictor, uThresh)
    init_state = human_traj{1};

    percents_over_time = [];
    for pt=1:length(preds)
        curr_pred = preds{pt};
        
        % Compute the likely enough states
        [opt_eps, P, X, Y] = ...
            compute_likely_states(curr_pred, predictor, uThresh);
        
        % Compute the FRS.
        frs_rad = v*dt*pt;
        
        % note in the case where pt == 1, it should just be no conservt.???
        
        states_occupied_in_frs = 0;
        total_states_in_frs = 0;
        for r = 1:predictor.rows
            for c = 1:predictor.cols
                d = sqrt((X(r,c) - init_state(1))^2 + (Y(r,c) - init_state(2))^2); 
                if d <= frs_rad
                    if P(r,c) == 1
                        states_occupied_in_frs = states_occupied_in_frs + 1;
                    end
                    total_states_in_frs = total_states_in_frs + 1;
                end
            end
        end
        percents_over_time(end+1) = states_occupied_in_frs/total_states_in_frs;
    end
    
    p_con = mean(percents_over_time);
end

%% Returns % of timesteps for which the true human state is contained 
%  within our predictions. The higher the better.
function p_acc = accuracy(preds, xcurr, tcurr, human_traj, predictor, uThresh)

    num_times_human_state_in_preds = 0;
    traj_idx = tcurr;
    for pt=1:length(preds)
        curr_pred = preds{pt};
        
        % Compute the likely enough states
        [~, P, X, Y] = ...
            compute_likely_states(curr_pred, predictor, uThresh);
        
        % Get the true state where the person actually went at this time.
        true_human_state = human_traj{traj_idx};
        
        % Compute the closest index in the prediction grid to this state.
        [i, j] = predictor.realToSim(true_human_state);
        
        % If human was predicted here, count it!
        if P(i,j) > 0
            num_times_human_state_in_preds = ...
                num_times_human_state_in_preds + 1;
        end
        traj_idx = traj_idx + 1;
    end
    
    p_acc = num_times_human_state_in_preds/length(preds);
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

