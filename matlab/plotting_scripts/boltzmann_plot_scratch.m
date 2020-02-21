clc
clear all

%% Grid
grid_min = [-2; -2; -0.1];  % Lower corner of computation domain
grid_max = [2; 2; 1.1];     % Upper corner of computation domain
N = [41; 41; 41];           % Number of grid points per dimension
g = createGrid(grid_min, grid_max, N);

%% Compute P(u_t | x_t, \beta = 0) = 
% \frac{e^{-v \cdot \Delta t - \|x_t + \Delta t \cdot f(x_t, u_t) - \theta\|_2}}{\sum_{u \in \mathcal{U}} e^{-v \cdot \Delta t - \|x_t + \Delta t \cdot f(x_t, u) - \theta\|_2}}

ind = 1;

figure;

%% Setup variables
x = g.xs; 

lb = -pi;
ub = pi;
v = 0.6;
numCtrls = 20;
inc = (ub-lb)/(numCtrls-1);
deltaT = 1;
theta = [0 0];

num_plots = 5;

ctrls = linspace(lb, ub, num_plots);

% Loop through controls
for ctrl_index = 1:num_plots
    
    u = ctrls(ctrl_index);
    
    %% Compute integral: \sum_{u \in \mathcal{U}} e^{-v \cdot \Delta t - \|x_t + \Delta t \cdot f(x_t, u) - \theta\|_2}
    % Note: Ignoring -v \cdot \Delta t as it is constant cancels out in expression
    
    intControls = 0.0;
    for i = 1:numCtrls

        % Find control
        u_temp = lb + (i-1)*inc;

        % Find next x by forward euler
        x1 = x{1} + deltaT .* v .* cos(u_temp);
        x2 = x{2} + deltaT .* v .* sin(u_temp);

        % Evaluate distance of next x to goal theta under L2 norm
        nrm = ((x1 - theta(1)).^2 + (x2 - theta(2)).^2).^(0.5);

        % Calculate value in summation: e^{-\| (x_t + \Delat t f(x_t,u_t)) - \theta \|_2}
        val = exp(1).^(-1.*nrm);

        % Add to running value of summation
        intControls = intControls + inc*val;
    end
    
    %% Computer e^{-v \cdot \Delta t - \|x_t + \Delta t \cdot f(x_t, u_t) - \theta\|_2}

    % x_{t+1} = x_t + \Delta t * f(x_t)
    x1 = x{1} + deltaT .* v .* cos(u);
    x2 = x{2} + deltaT .* v .* sin(u);

    nrm = ((x1 - theta(1)).^2 + (x2 - theta(2)).^2).^(0.5);
    val = exp(1).^(-1.*nrm);

    prob_u = val ./ intControls;
    
    % Are the x-y directions incorrect
    prob_u_high = prob_u(:,:,38); % high posterior prob of ~1.0
    
    %% Plot
    subplot(2,3,ind)
    heatmap(prob_u_high);
    title(["u = " u])
    
    ind = ind + 1;
end
