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

lb = -pi+1e-2;
ub = pi;
v = 0.6;
numCtrls = 100;
inc = (ub-lb)/(numCtrls-1);
deltaT = 1;
theta = [1 1];

num_plots = 5;

ctrls = linspace(lb, ub, num_plots);

% Get a 2D grid structure for just x and y values.
grid2D = createGrid(grid_min(1:2), grid_max(1:2), N(1:2));

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
        val = exp(-1.*nrm); %exp(1).^(-1.*nrm);

        % Add to running value of summation
        intControls = intControls + inc*val;
    end
    
    %% Computer e^{-v \cdot \Delta t - \|x_t + \Delta t \cdot f(x_t, u_t) - \theta\|_2}

    % x_{t+1} = x_t + \Delta t * f(x_t)
    x1 = x{1} + deltaT .* v .* cos(u);
    x2 = x{2} + deltaT .* v .* sin(u);

    nrm = ((x1 - theta(1)).^2 + (x2 - theta(2)).^2).^(0.5);
    val = exp(-1.*nrm); %exp(1).^(-1.*nrm);

    prob_u = val ./ intControls;
    
    % NOTE: The prob_u is actually repeated for each 3rd index 
    %   so we can take any of the entries for plotting purposes
    %  (to see this check out a 5x5x5 grid and look at what prob_u looks
    %  like. This is beacuse its not actually corresponding to different
    %  posteriors, since we are only trying to plot our observation model's
    %  results, P(u | x, theta) and not the posterior P(theta | u, x). 
    
    % Should this be transposed to corresponed to x-y grid? If so, then
    % graphs make sense
    % prob_u_high = prob_u(:,:,end); % high posterior prob of ~1.0
    
    %% Plot
    subplot(2,3,ind)
    
    % Generate pseudocolor map of probabilities of u at each state in grid.
    h = pcolor(grid2D.xs{1}, grid2D.xs{2}, prob_u(:,:,1));
    colorbar
    caxis([0,0.5]);
    set(h, 'EdgeColor', 'none');
    title(["u = " u])
    
    hold on 
    % Plot the goal location.
    scatter([theta(1)], [theta(2)], 25, 'r', 'filled');
    
    ind = ind + 1;
end