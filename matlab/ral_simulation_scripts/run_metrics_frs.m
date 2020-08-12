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

% Grid cell size.
r_real = 0.1;

% Velocity
v = 0.6;

% Time info.
dt = r_real / v;
T = 2;                          % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)
discrete_times = 0:1:H+1;       % (used for data storing and plotting)

%% Data tracking.
all_con = [];
all_acc = [];

sim_cond_idx = 1;

% timing. 
one_sim_start_t = tic;

avg_con_per_sigma = [];
avg_acc_per_sigma = [];

%% Run metrics!
for sigma = human_sigmas
    traj_map = all_trajs(num2str(sigma));
    this_sigma_con = [];
    this_sigma_acc = [];
    for init_idx = 1:1 %init_state = human_init_conds
        %x0 = init_state{1};
        
        x0 = human_init_conds{init_idx};
        human_traj = traj_map(num2str(x0));
        
        % Main simulation loop!
        acc_rn = [];
        con_rn = [];
        for t=1:simT
            xcurr = human_traj{t}{1};
            ucurr = human_traj{t}{2};
            
            %% Predict the human.
            fprintf('FRS %d / %d at sim t = %d / %d...\n', ...
                sim_cond_idx, length(human_sigmas), ...
                t, simT);
            
            fprintf('Computing metrics...\n');
  
            %% Compute % overconservative. 
            cons_debug = false;
            p_con = 1.0; %conservativism(xcurr, v, dt, H);
            
            %% Compute % accurate. 
            acc_debug = false;
            p_acc = accuracy(xcurr, t, human_traj, v, dt, H);
            
            %% Store metrics.
            all_con(end+1) = p_con;
            all_acc(end+1) = p_acc;
            
            % Store how predictions are doing at each simulation setp.
            acc_rn(end+1) = p_acc;
            con_rn(end+1) = p_con;
            
        end
        
        sim_cond_idx = sim_cond_idx + 1;
        
        % Average over simulation for this intiial condition and sigma. 
        this_sigma_con(end+1) = mean(con_rn);
        this_sigma_acc(end+1) = mean(acc_rn);
        
        fprintf('for sigma=%f, init=(%f,%f): \n', sigma, x0(1), x0(2));
        fprintf('    AVG CON = %f\n', mean(con_rn));
        fprintf('    AVG ACC = %f\n', mean(acc_rn));
    end
    
    % Average over all initial conditions for this sigma. 
    avg_con_per_sigma(end+1) = mean(this_sigma_con);
    avg_acc_per_sigma(end+1) = mean(this_sigma_acc);

end

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
% function p_con = conservativism(preds, xcurr, v, dt, ...
%     predictor, uThresh, debug)
%     
%     if debug
%         figure
%     end
% 
%     init_state = xcurr; % human_traj{1};
% 
%     percents_over_time = [];
%     for pt=1:length(preds)
%         curr_pred = preds{pt};
%         
%         % Compute the likely enough states
%         [~, P, X, Y] = ...
%             compute_likely_states(curr_pred, predictor, uThresh);
%         
%         % Compute the FRS.
%         frs_rad = v*dt*(pt-1);
%         
%         if debug
%             all_plt_handles = {};
%         end
% 
%         states_occupied_in_frs = 0;
%         total_states_in_frs = 0;
%         for r = 1:predictor.rows
%             for c = 1:predictor.cols
%                 d = sqrt((X(r,c) - init_state(1))^2 + (Y(r,c) - init_state(2))^2); 
%                 % In case of FRS at timestep = 0, the only valid state is
%                 % the intial state of the human.               
%                 if d <= frs_rad
%                     if debug
%                         hold on
%                         % Plot prediction contour.
%                         if P(r,c) == 1
%                             sp = scatter(X(r,c), Y(r,c), 20, 'b', 'filled');
%                         else
%                             sp = scatter(X(r,c), Y(r,c), 20, 'r');
%                         end
%                         % Plot state of human.
%                         s = scatter(init_state(1), init_state(2), 10, 'k', 'filled');
% 
%                         pos = [init_state(1)-frs_rad, init_state(2)-frs_rad, frs_rad*2, frs_rad*2];
%                         rec = rectangle('Position',pos,'Curvature',1);
%                         xlim([-4, 4]);
%                         ylim([-4, 4]);
%                         set(gcf,'Position',[100 100 700 700]);
%                         set(gcf,'color','w');
%                         whitebg('w');
%                         all_plt_handles{end+1} = s;
%                         all_plt_handles{end+1} = sp;
%                         all_plt_handles{end+1} = rec;
%                         grid on
%                         hold off
%                     end
%                         
%                     if P(r,c) == 1
%                         states_occupied_in_frs = states_occupied_in_frs + 1;
%                     end
%                     total_states_in_frs = total_states_in_frs + 1;
%                 end
%             end
%         end
%         
%         if debug 
%             hold on
%             [M, c] = contour(X, Y, P, [1, 1]);
%             c.LineWidth = 2;
%             c.EdgeColor = [0,0,1];
%             hold off
%             for hand = all_plt_handles
%                 delete(hand{1});
%             end
%         end
%         
%         if total_states_in_frs == 0 % case where frs_rad == 0
%             total_states_in_frs = 1;
%         end
%         percents_over_time(end+1) = states_occupied_in_frs/total_states_in_frs;
%     end
%     
%     p_con = mean(percents_over_time);
% end

%% Returns % of timesteps for which the true human state is contained 
%  within our predictions. The higher the better.
function p_acc = accuracy(xcurr, tcurr, human_traj, v, dt, H)

    num_times_human_state_in_preds = 0;
    traj_idx = tcurr;
    for pt=1:H
        
         % Compute the FRS radius.
        frs_rad = v*dt*(pt-1);
        
        % Get the true state where the person actually went at this time.
        true_human_state = human_traj{traj_idx}{1};
        
        % Compute distance of the true human state at this time from the
        % initial state. 
        d = sqrt((true_human_state(1) - xcurr(1))^2 + (true_human_state(2) - xcurr(2))^2); 
        
        % If human was predicted here, count it!
        if d <= frs_rad + 1e10
            num_times_human_state_in_preds = ...
                num_times_human_state_in_preds + 1;
        else
            % Plot state of human.
            figure
            hold on
            s = scatter(true_human_state(1), true_human_state(2), 40, 'k', 'filled');

            cs = scatter(xcurr(1), xcurr(2), 10, 'k');
            pos = [xcurr(1)-frs_rad, xcurr(2)-frs_rad, frs_rad*2, frs_rad*2];
            rec = rectangle('Position',pos,'Curvature',1);
            xlim([-4, 4]);
            ylim([-4, 4]);
            hold off
        end
        traj_idx = traj_idx + 1;
        traj_idx = min(length(human_traj), traj_idx);
    end
    
    % Option 1:
    % Return ratio of accurate predictions to pred horizon.
    % p_acc = num_times_human_state_in_preds/length(preds);
    
    % Option 2: 
    % Return 1 if all true human states predicted, 0 otherwise.
    if num_times_human_state_in_preds == H
        p_acc = 1;
    else
        p_acc = 0;
    end
end


