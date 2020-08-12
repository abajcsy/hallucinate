clear all
close all

%% Random seed.
rng(13);

%% Simulated human params.

% Velocity
v = 0.6;

% Control bounds
uRange = [-pi; pi];

% Known human goal locations. 
goals = {[2; 2], [2; -2]}; 

% True goal that human is going to.
trueGoal = goals(1);

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

% Grid cell size.
r = 0.1;
dt = r / v;

%% Plot for debugging?
plot = false;

%% Number of cells to divide x and y into.
num_init_cells = 4;
init_res_x = (grid_max(1) - grid_min(1))/num_init_cells;
init_res_y = (grid_max(2) - grid_min(2))/num_init_cells;

x = grid_min(1):init_res_x:grid_max(1);
y = grid_min(2):init_res_y:grid_max(2);
cell_bounds = {};
for yi = y(1:end-1)
    for xi = x(1:end-1)
        if (xi == goals{1}(1) && yi == goals{1}(2)) || ...
                xi == goals{2}(1) && yi == goals{2}(2)
            continue;
        end
        cell_bounds{end+1} = [xi; yi; xi+init_res_x; yi+init_res_x];
    end
end

%% Setup simulation conditions.

% Set of standard deviations of the simulated human action.
human_sigmas = [0, pi/4, pi/2, 3*pi/4, pi];

% (Randomly sample) set of initial conditions from where to simulate the human. 
human_init_conds = {};
for i=1:length(cell_bounds)
    bounds = cell_bounds{i};
    % generate random number in the interval (a,b) with 
    % the formula r = a + (b-a) .* rand(1,1).
    rx = bounds(1) + (bounds(3)-bounds(1)) * rand(1,1);
    ry = bounds(2) + (bounds(4)-bounds(2)) * rand(1,1);
    human_init_conds{end+1} = [rx; ry];
    
    if plot
        hold on;
        scatter(rx, ry, 'filled');
    end
end

if plot
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    set(gca,'xtick',x);
    set(gca,'ytick',y);
    grid on
end

%% Setup simulation.

% Number of simulation steps. (real sim time = T*dt)
simT = 102; % (i.e. 17 seconds)

% Create map to store for each sigma, for each initial condition the
% trajectory the simulated human took:
%       key: sigma, value: Map --> key: init_cond, value: trajectory
%                                                         [state, action]
all_trajs = containers.Map;

%% Generate the simulated data!
for sigma = human_sigmas
    traj_map = containers.Map;
    for init_state = human_init_conds
        x0 = init_state{1};
        
        % Create simulated human.
        human = GaussianGoalHuman(x0, v, sigma, uRange, trueGoal);
        
        % State. 
        xcurr = x0;
        human_states = {};
        times = [];
        
        if plot
            figure
            hold on
        end
        
        % Main simulation loop!
        for t=1:simT
            if plot
                scatter(xcurr(1), xcurr(2), 'filled');
            end 
           
            
            % Forward simulate human. 
            [xnext, ucurr] = human.simulate(xcurr, t*dt, dt);
            while xnext(1) < grid_min(1) || xnext(1) > grid_max(1) || ...
                    xnext(2) < grid_min(2) || xnext(2) > grid_max(2)
                % If the next state is outside the grid, resample.
                [xnext, ucurr] = human.simulate(xcurr, t*dt, dt);
            end
            
            % Store current state.
            human_states{end+1} = {xcurr, ucurr};
            times(end+1) = (t-1)*dt;
            
            % Update state.
            xcurr = xnext;
        end
        
        % Save this trajectory.
        traj_map(num2str(x0)) = human_states;
        
        if plot
            scatter(trueGoal{1}(1), trueGoal{1}(2), 'r');
            xlim([grid_min(1), grid_max(1)]);
            ylim([grid_min(2), grid_max(2)]);
            set(gca,'xtick',x);
            set(gca,'ytick',y);
            grid on
            pause(0.5);
        end
    end
    % Save the trajectories for all initial conditions for this sigma.
    all_trajs(num2str(sigma)) = traj_map;
end
repo = what('hallucinate');
filename = strcat('sigma_init_cond.mat');
save(strcat(repo.path, '/ral_data_revise/human_traj_data/', filename), ...
    'all_trajs', 'human_sigmas', 'human_init_conds', 'simT', 'dt');
