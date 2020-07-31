clear all
close all

% TODO
% -- make sure that the simulated human can't go outside of the compute
% grid
% -- add in saving functionality for data

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

% Plot for debugging?
plot = true;

% Number of cells to divide x and y into.
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
simT = 80;

%% Generate the simulated data!
for sigma = human_sigmas
    bla = 1;
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
            
            % Store current state.
            human_states{end+1} = xcurr;
            times(end+1) = (t-1)*dt;
            
            % Forward simulate human. 
            [xnext, ucurr] = human.simulate(xcurr, t*dt, dt);
            xcurr = xnext;
        end
        
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
end
