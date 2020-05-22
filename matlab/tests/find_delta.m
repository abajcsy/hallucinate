clear all
close all

grid_min = [-8; -8; -0.1];  % Lower corner of computation domain
grid_max = [8; 8; 1.1];     % Upper corner of computation domain
N = [51; 51; 51];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);
xs = g.xs;

% control bounds.
urange = [-pi, pi];

% goal locations (meters).
g1 = [1; 0];

v = 0.2;
deltaT = 1;

% number of discrete controls.
num_ctrls = 51;

% get discretized controls.
us = gen_controls(urange, num_ctrls);

% action probabilities
us_probs_0 = cell(1,length(us));
% us_probs_0_1 = zeros(1,length(us));
us_probs_1 = cell(1,length(us));
us_probs = cell(1,length(us));

prior = xs{3};

% get the probability of each action.
for i=1:num_ctrls
    u = us(i);
    
    % Compute action probabilities
    pu_0 = pu_beta(u,xs,0,us,v,deltaT,g1);
    pu_1 = pu_beta(u,xs,1,us,v,deltaT,g1);
    pu = pu_0 .* (1-prior) + pu_1 .* prior;
    
    % plotting info...
    us_probs_0{i} = pu_0;
    us_probs_1{i} = pu_1;
    us_probs{i} = pu;
end

range = find_prob_range(us_probs_1);
range

%% Find probability range of controls
function range = find_prob_range(us_probs)
    probs = {};
    for i=1:length(us_probs)
        u_probs = us_probs(i);

        % Convert into an N1 x N2 x N3 x numCtrls array
        if i == 1
            probs = u_probs;
        else
            probs = cat(4, probs, u_probs);
        end
    end
    prob_matrix = cell2mat(probs);
    
    min_prob = min(prob_matrix,[],4);
    max_prob = max(prob_matrix,[],4);
    
    delta_max = max(max_prob,[],'all');
    delta_min = max(min_prob,[],'all');
    
    range = [delta_min; delta_max];
end

%% Generate discrete num_ctrls ranging from [urange(1), urange(2)) 
function us = gen_controls(urange, num_ctrls)
    incr = (urange(2) - urange(1))/num_ctrls;
    us = zeros(1,num_ctrls);
    u = urange(1);
    for i=1:num_ctrls
        us(i) = u;
        u = u + incr;
    end
end

%% Find probability of action for specific beta
function pu = pu_beta(u,x,beta,us,v,deltaT,goal)
    us_shape = size(us);
    % Calculta sum of 
    sumControls = 0.0;
    for i = 1:us_shape(2)

        % Get discrete control.
        u_i = us(i);

        % Compute the Q-value of each state and control.
        qval = qFunction(x, u_i, goal, v, deltaT);

        % Calculate value in summation: 
        %   e^{-||(x_t + \Deltat t f(x_t,u_t)) - \theta||_2}
        u_i_val = exp(beta .* qval);

        % Add to running value of summation
        sumControls = sumControls + u_i_val; 
    end
    u_val = exp(beta .* qFunction(x, u, goal, v, deltaT));
    pu = u_val ./ sumControls;
end

function qval = qFunction(x, u, theta, v, deltaT)
    % Find next x by forward euler
    x1 = x{1} + deltaT * v * cos(u);
    x2 = x{2} + deltaT * v * sin(u);

    % Evaluate distance of next x to goal theta under L2 norm
    qval = -1*((x1 - theta(1)).^2 + (x2 - theta(2)).^2).^(1.);
end