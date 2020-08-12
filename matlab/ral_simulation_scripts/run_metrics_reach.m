clear all
close all

%% Load data.
repo = what('hallucinate');
filename = 'sigma_init_cond.mat';
load(strcat(repo.path, '/ral_data_revise/human_traj_data/', filename));

%% (BA-FRS) Predictor params.

% Grid
grid_min = [-4; -4; -0.1];  % Lower corner of computation domain
grid_max = [4; 4; 1.1];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

% gamma in continuous-time P(beta = 0) dynamics
gamma = 0.4; %1;

% Number of discrete controls
numCtrls = 31;

% Are we using dynamic of static beta model?
betaModel = 'static';
% (No) Dynamic beta parameters
extraArgs = [];

% (reachablity) Threshold to determine likely controls -- used for
% plotting.
uThresh = 0.04; %0.03; % 0.02;

% Variance on Gaussian observation model.
sigma = pi/4;

% Prior. 
prior = [0.5, 0.5];
Pgoal1 = prior(1); 

% Velocity
v_real = 0.6;
v_pred = v_real; %0.8638;

% Control bounds
uRange = [-pi; pi];

% Known human goal locations. 
goals = {[2; 2], [2; -2]}; 

% Setup dynamical system (i.e. predictor)
z0 = [0; 0; Pgoal1]; %NOTE: this is incorrect and is overriden in sim. 
predictor = GaussianTwoGoalHuman(z0, v_pred, uRange, gamma, goals, ...
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
tMax = 2; 
dt = 0.1667;
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
% HJIextraArgs.visualize.valueSet = 1;
% HJIextraArgs.visualize.initialValueSet = 0;
% HJIextraArgs.visualize.figNum = 1; 
% HJIextraArgs.visualize.deleteLastPlot = false; 
% HJIextraArgs.visualize.viewGrid = true;
% HJIextraArgs.visualize.viewAxis = [grid_min(1) grid_max(1) ...
%                                    grid_min(2) grid_max(2) ...
%                                    -0.1 1.1];
% HJIextraArgs.visualize.xTitle = '$p^x$';
% HJIextraArgs.visualize.yTitle = '$p^y$';
% HJIextraArgs.visualize.zTitle = '$P(goal_1)$';
% HJIextraArgs.visualize.fontSize = 15;

HJIextraArgs.ignoreBoundary = 0; 
HJIextraArgs.quiet = 1; 

% Compute FRS.
minWith = 'zero';

%% Human footprint. 
human_footprint_rad = 0.3/2; % radius of footprint.

%% End setup.

all_con = [];
all_acc = [];

sim_cond_idx = 1;

% timing. 
one_sim_start_t = tic;

avg_con_per_sigma = [];
avg_acc_per_sigma = [];

%% Pick which sigma and init condition to debug.
% sigma = human_sigmas(1);
% x0 = human_init_conds{1};
% 
% traj_map = all_trajs(num2str(sigma));
% human_traj = traj_map(num2str(x0));


%% Run metrics!
for sigma = human_sigmas
    traj_map = all_trajs(num2str(sigma));
    this_sigma_con = [];
    this_sigma_acc = [];
    for init_idx = 1:1 %init_state = human_init_conds
        
        x0 = human_init_conds{init_idx};
        human_traj = traj_map(num2str(x0));
        
        % Main simulation loop!
        acc_rn = [];
        con_rn = [];
        for t=1:simT
            xcurr = human_traj{t}{1};
            ucurr = human_traj{t}{2};

            %% Update initial set and predictor info.
            z0 = [xcurr(1); xcurr(2); Pgoal1];
            data0 = shapeSphere(g, z0, R);
            predictor.x = z0;
            schemeData.dynSys = predictor;

            %% Predict!
            fprintf('Predicting %d / %d at sim t = %d / %d...\n', ...
                sim_cond_idx, length(human_sigmas), ...
                t, simT);
            [preds, pred_times, ~] = ...
                HJIPDE_solve_pred(data0, tau, schemeData, minWith, HJIextraArgs);

            fprintf('Computing metrics...\n');

            % Compute % overconservative. 
            cons_debug = false;
            p_con = conservativism(preds, xcurr, v_real, dt, g, cons_debug);

            % Compute % accurate. 
            acc_debug = false;
            p_acc = accuracy(preds, xcurr, t, human_traj, g, human_footprint_rad, v_real, dt, acc_debug);

            all_con(end+1) = p_con;
            all_acc(end+1) = p_acc;

            % Store how predictions are doing at each simulation setp.
            acc_rn(end+1) = p_acc;
            con_rn(end+1) = p_con;

            %% Update Posterior.
            xcurr_cell = {xcurr(1), xcurr(2)}; 
            pu_goal1 = predictor.PuGivenGoal_normalized(ucurr, xcurr_cell, 1);
            pu_goal2 = predictor.PuGivenGoal_normalized(ucurr, xcurr_cell, 2);
            Pgoal1 = (pu_goal1*Pgoal1)/(pu_goal1*Pgoal1 + pu_goal2*(1-Pgoal1));
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

fprintf('======= [BA-FRS] TOTAL ELAPSED TIME = %f ======= \n', one_sim_end_t);
fprintf('======= [BA-FRS] TOTAL AVG CONSERVATIVE = %f ======= \n', avg_con);
fprintf('======= [BA-FRS] TOTAL AVG ACCURATE = %f ======= \n', avg_acc);
for si = 1:length(human_sigmas)
    fprintf('       FOR sigma = %f, AVG CONSERVATIVE = %f ======= \n', ...
        human_sigmas(si), avg_con_per_sigma(si));
    fprintf('       FOR sigma = %f, AVG ACCURACY = %f ======= \n', ...
        human_sigmas(si), avg_acc_per_sigma(si));
    fprintf('------------------------------------------------------\n');
end

%% Returns percent of the full FRS that the predictions occupy.
% The lower the better.
function p_con = conservativism(preds, xcurr, v, dt, g, debug)
    
    if debug
        figure
    end

    init_state = xcurr; 

    percents_over_time = [];
    dimsToRemove = [0 0 1];
    num_preds = length(preds(1,1,1,:));
    for pt=1:num_preds
        
        % Project over all beliefs.
        [g2D, preds2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
                
        % Compute the FRS.
        frs_rad = v*dt*(pt-1);
        
        % Get distance 
        d = sqrt((g2D.xs{1} - init_state(1)).^2 + (g2D.xs{2} - init_state(2)).^2);
        
        idxs_in_frs = find(d <= frs_rad);
        total_states_in_frs = length(idxs_in_frs);
        num_pred_states_in_frs = 1*(preds2D(idxs_in_frs) <= 0);
        states_occupied_in_frs = sum(num_pred_states_in_frs);
       
        if total_states_in_frs == 0 % case where frs_rad == 0
            total_states_in_frs = 1;
        end
        percents_over_time(end+1) = states_occupied_in_frs/total_states_in_frs;
    end
    
    p_con = mean(percents_over_time);
end

%% Returns % of timesteps for which the true human state is contained 
%  within our predictions. The higher the better.
function p_acc = accuracy(preds, xcurr, tcurr, human_traj, g, human_footprint_rad, v, dt, debug)

    if debug
        figure
    end

    num_times_human_state_in_preds = 0;
    traj_idx = tcurr;
    
    dimsToRemove = [0 0 1];
    num_preds = length(preds(1,1,1,:));
    for pt=1:num_preds
        % Project over all beliefs.
        [g2D, preds2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
        
        % Get the true state where the person actually went at this time.
        true_human_state = human_traj{traj_idx}{1};
                
        % Get distance between true human state and grid.
        d_center = sqrt((g2D.xs{1} - true_human_state(1)).^2 + (g2D.xs{2} - true_human_state(2)).^2);
        d_human = d_center - human_footprint_rad;
        

        human_idxs = find(d_human <= 0);
        if human_footprint_rad == 0
            % Find the closest index in the prediction grid to this state.
            mval = min(d_human, [],'all');
            human_idxs = find(d_human == mval);
        end
        
        boundary_tol = 0; %0.01;
        
        % If human was predicted here, count it!
        if any(preds2D(human_idxs) <= boundary_tol) %preds2D(midx) <= boundary_tol
            num_times_human_state_in_preds = ...
                num_times_human_state_in_preds + 1;
        end
        traj_idx = traj_idx + 1;
        traj_idx = min(length(human_traj), traj_idx);
        
        if debug
            hold on
            visSetIm(g2D, preds2D);
            % real human state
            scatter(true_human_state(1), true_human_state(2), 'k', 'filled');
            % gridded human state
            scatter(g2D.xs{1}(human_idxs), g2D.xs{2}(human_idxs), 'r');
            xlim([-4, 4]);
            ylim([-4, 4]);
            hold off
            
%             hold on
%             visSetIm(g, preds(:,:,:,pt));
%             % Compute the FRS
%             R = 0; %0.1;
%             frs_rad = v*dt*(pt-1) + R;
%             pos = [ xcurr(1)-frs_rad,  xcurr(2)-frs_rad, frs_rad*2, frs_rad*2];
%             rec = rectangle('Position',pos,'Curvature',1);
%             scatter(2, 2, 40, 'r', 'filled');
%             scatter(2, -2, 40, 'r', 'filled');
%             xlim([-4, 4]);
%             ylim([-4, 4]);
%             zlim([0,1.1]);
%             grid on
%             hold off
        end
    end
    
    % Option 1:
    % Return ratio of accurate predictions to pred horizon.
    % p_acc = num_times_human_state_in_preds/length(preds);
    
    % Option 2: 
    % Return 1 if all true human states predicted, 0 otherwise.
    if num_times_human_state_in_preds == num_preds
        p_acc = 1;
    else
        p_acc = 0;
    end
end

