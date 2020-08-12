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
r_real = 0.1;
r_pred = 0.144;

% Velocity
v = 0.6;

% Known human goal locations. 
goals = {[2; 2], [2; -2]}; 

% Create the predictor.
predictor = LatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
    grid_min, grid_max, r_pred);

dt = r_real / v;
T = 2;                          % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)
discrete_times = 0:1:H+1;       % (used for data storing and plotting)

all_con = [];
all_acc = [];

sim_cond_idx = 1;

% timing. 
one_sim_start_t = tic;

avg_con_per_sigma = [];
avg_acc_per_sigma = [];

%% Human footprint. 
human_footprint_rad = 0.3/2; % radius of footprint.

%% Pick which sigma and init condition to debug.
sigma = human_sigmas(end);
x0 = human_init_conds{1};

%% Run metrics!

traj_map = all_trajs(num2str(sigma));
human_traj = traj_map(num2str(x0));

this_sigma_con = [];
this_sigma_acc = [];

%% Main simulation loop!
acc_rn = [];
con_rn = [];
for t=1:simT
    xcurr = human_traj{t}{1};
    ucurr = human_traj{t}{2};

    % Predict the human.
%             fprintf('Predicting %d / %d at sim t = %d / %d...\n', ...
%                 sim_cond_idx, length(human_sigmas)*length(human_init_conds), ...
%                 t, simT);
    fprintf('Predicting %d / %d at sim t = %d / %d...\n', ...
        sim_cond_idx, length(human_sigmas), ...
        t, simT);
    preds = predictor.predict(xcurr, H);

    fprintf('Computing metrics...\n');

    % Compute % overconservative. 
    cons_debug = false;
    p_con = conservativism(preds, xcurr, human_traj, v, dt, ...
                            predictor, uThresh, cons_debug);

    % Compute % accurate. 
    acc_debug = true;
    p_acc = accuracy(preds, xcurr, t, human_traj, ...
                      predictor, uThresh, human_footprint_rad, acc_debug);

    %fprintf('for sigma=%f, init=(%f,%f): \n', sigma, x0(1), x0(2));
%             fprintf('    CON = %f\n', p_con);
%             fprintf('    ACC = %f\n', p_acc);

    all_con(end+1) = p_con;
    all_acc(end+1) = p_acc;

    % Store how predictions are doing at each simulation setp.
    acc_rn(end+1) = p_acc;
    con_rn(end+1) = p_con;
    
    %% Update Posterior.
    [i,j] = predictor.realToSim(xcurr);
    s = [i,j];
    pu_goal1 = predictor.Pu_given_x_g_normalized(ucurr, s, 1);
    pu_goal2 = predictor.Pu_given_x_g_normalized(ucurr, s, 2);
    Pgoal1 = (pu_goal1*Pgoal1)/(pu_goal1*Pgoal1 + pu_goal2*(1-Pgoal1));

    % Update predictor info.
    new_prior = [Pgoal1, (1-Pgoal1)];
    predictor.prior = containers.Map([1:length(new_prior)], new_prior);
end

sim_cond_idx = sim_cond_idx + 1;

% Average over simulation for this intiial condition and sigma. 
this_sigma_con(end+1) = mean(con_rn);
this_sigma_acc(end+1) = mean(acc_rn);

fprintf('for sigma=%f, init=(%f,%f): \n', sigma, x0(1), x0(2));
fprintf('    AVG CON = %f\n', mean(con_rn));
fprintf('    AVG ACC = %f\n', mean(acc_rn));

% Average over all initial conditions for this sigma. 
avg_con_per_sigma(end+1) = mean(this_sigma_con);
avg_acc_per_sigma(end+1) = mean(this_sigma_acc);


one_sim_end_t = toc(one_sim_start_t);

% Return the average over all accuracy and conservativism!
avg_con = mean(all_con);
avg_acc = mean(all_acc);

fprintf('======= TOTAL ELAPSED TIME = %f ======= \n', one_sim_end_t);
fprintf('======= TOTAL AVG CONSERVATIVE = %f ======= \n', avg_con);
fprintf('======= TOTAL AVG ACCURATE = %f ======= \n', avg_acc);
for si = 1:length(human_sigmas)
    fprintf('       FOR sigma = %f, AVG CONSERVATIVE = %f ======= \n', ...
        human_sigmas(si), avg_con_per_sigma(si));
    fprintf('       FOR sigma = %f, AVG ACCURACY = %f ======= \n', ...
        human_sigmas(si), avg_acc_per_sigma(si));
    fprintf('------------------------------------------------------\n');
end

%% Returns percent of the full FRS that the predictions occupy.
% The lower the better.
function p_con = conservativism(preds, xcurr, human_traj, v, dt, ...
    predictor, uThresh, debug)
    
    if debug
        figure
    end

    init_state = xcurr; % human_traj{1};

    percents_over_time = [];
    for pt=1:length(preds)
        curr_pred = preds{pt};
        
        % Compute the likely enough states
        [~, P, X, Y] = ...
            compute_likely_states(curr_pred, predictor, uThresh);
        
        % Compute the FRS.
        frs_rad = v*dt*(pt-1);
        
        if debug
            all_plt_handles = {};
        end

        states_occupied_in_frs = 0;
        total_states_in_frs = 0;
        for r = 1:predictor.rows
            for c = 1:predictor.cols
                d = sqrt((X(r,c) - init_state(1))^2 + (Y(r,c) - init_state(2))^2); 
                % In case of FRS at timestep = 0, the only valid state is
                % the intial state of the human.               
                if d <= frs_rad
                    if debug
                        hold on
                        % Plot prediction contour.
                        if P(r,c) == 1
                            sp = scatter(X(r,c), Y(r,c), 20, 'b', 'filled');
                        else
                            sp = scatter(X(r,c), Y(r,c), 20, 'r');
                        end
                        % Plot state of human.
                        s = scatter(init_state(1), init_state(2), 10, 'k', 'filled');

                        pos = [init_state(1)-frs_rad, init_state(2)-frs_rad, frs_rad*2, frs_rad*2];
                        rec = rectangle('Position',pos,'Curvature',1);
                        xlim([-4, 4]);
                        ylim([-4, 4]);
                        set(gcf,'Position',[100 100 700 700]);
                        set(gcf,'color','w');
                        whitebg('w');
                        all_plt_handles{end+1} = s;
                        all_plt_handles{end+1} = sp;
                        all_plt_handles{end+1} = rec;
                        grid on
                        hold off
                    end
                        
                    if P(r,c) == 1
                        states_occupied_in_frs = states_occupied_in_frs + 1;
                    end
                    total_states_in_frs = total_states_in_frs + 1;
                end
            end
        end
        
        if debug 
            hold on
            [M, c] = contour(X, Y, P, [1, 1]);
            c.LineWidth = 2;
            c.EdgeColor = [0,0,1];
            hold off
            for hand = all_plt_handles
                delete(hand{1});
            end
        end
        
        if total_states_in_frs == 0 % case where frs_rad == 0
            total_states_in_frs = 1;
        end
        percents_over_time(end+1) = states_occupied_in_frs/total_states_in_frs;
    end
    
    p_con = mean(percents_over_time);
end

%% Returns % of timesteps for which the true human state is contained 
%  within our predictions. The higher the better.
function p_acc = accuracy(preds, xcurr, tcurr, human_traj, ...
        predictor, uThresh, human_footprint_rad, debug)

    if debug
        figure
    end

    num_times_human_state_in_preds = 0;
    traj_idx = tcurr;
    for pt=1:length(preds)
        curr_pred = preds{pt};
        
        % Compute the likely enough states
        [~, P, X, Y] = ...
            compute_likely_states(curr_pred, predictor, uThresh);
        
        % Get the true state where the person actually went at this time.
        true_human_state = human_traj{traj_idx}{1};
        
        % Compute the closest index in the prediction grid to this state.
        [i, j] = predictor.realToSim(true_human_state);
        
        % Get distance between true human state and grid.
        d_center = sqrt((X - true_human_state(1)).^2 + (Y - true_human_state(2)).^2);
        % Account for the human footprint.
        d_human = d_center - human_footprint_rad;
       
        human_idxs = find(d_human <= 0);
        xreal_foot = X(human_idxs);
        yreal_foot = Y(human_idxs);
        num_predicted = 0;
        for ns = 1:length(human_idxs)
            [ir,jr] = predictor.realToSim([xreal_foot(ns); yreal_foot(ns)]);
            if P(ir,jr) > 0
                num_predicted = num_predicted + 1;
            end
        end
        
        % If human was predicted here, count it!
        if num_predicted > 0
            num_times_human_state_in_preds = ...
                num_times_human_state_in_preds + 1;
        end
        
        traj_idx = traj_idx + 1;
        traj_idx = min(length(human_traj), traj_idx);
        
        if debug
            hold on
            % Plot prediction contour.
            color = [(length(preds)-pt)/length(preds),0,1];
            color_h = [(length(preds)-pt)/length(preds), ...
                        (length(preds)-pt)/length(preds), ...
                        (length(preds)-pt)/length(preds)];
            [M, c] = contour(X, Y, P, [1, 1]);
            c.LineWidth = 2;
            c.EdgeColor = color;
            % Plot state of human.
            s = scatter(true_human_state(1), true_human_state(2), 40, color_h, 'filled');
            s.MarkerEdgeColor = color;
            s.LineWidth = 1.5;
            
            % Plot GRIDDED state of human.
%             [xr, yr] = predictor.simToReal([i,j]);
%             sg = scatter(xr, yr, 40, 'r');
            sg = scatter(X(human_idxs), Y(human_idxs), 'r');
            xlim([-4, 4]);
            ylim([-4, 4]);
            set(gcf,'Position',[100 100 700 700]);
            set(gcf,'color','w');
            title(strcat("t0 = ", num2str(tcurr), ", k = ", num2str(pt)));
            whitebg('w');
            delete(c);
            delete(s);
            delete(sg);
            grid on
            hold off
        end
    end
    
    % Option 1:
    % Return ratio of accurate predictions to pred horizon.
    % p_acc = num_times_human_state_in_preds/length(preds);
    
    % Option 2: 
    % Return 1 if all true human states predicted, 0 otherwise.
    if num_times_human_state_in_preds == length(preds)
        p_acc = 1;
    else
        p_acc = 0;
    end
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

