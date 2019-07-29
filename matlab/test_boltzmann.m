function test_boltzmann()
clf
clear all

%% Setup.
% Define the action set.
global UP DOWN LEFT RIGHT
UP = 1;
RIGHT = 2;
DOWN = 3;
LEFT = 4;

% Initial (x,y) state of the human.
x0 = [0;0];

% Known goal.
g = [0;4];

% Set of values beta can take.
betas = [0, 0.1, 1, 100];

% Prior over beta.
P_b = [0.25,0.25,0.25,0.25];

%% Update prior over \beta given new measurement.
% Simulate a bunch of "good" actions (UP) and "bad" actions (DOWN).
P_b_new = P_b;

num_good_us = 100;
num_bad_us = 100;

for i=1:num_good_us
    P_b_new = update_Pb(UP, P_b_new, x0, g, betas);
end

for i=1:num_bad_us
    P_b_new = update_Pb(DOWN, P_b_new, x0, g, betas);
end

%% Predict next state distribution given your posterior.

% predict forward 1 timestep.
P_x1 = update_Px(x0, g, betas, P_b_new);

%% Plotting. 
% Plot the posterior over beta.
figure(1)
heatmap(P_b_new, 'ColorLimits',[0 1]);
ax = gca;
ax.XData = ["b = 0" "0.1" "1" "10"];
ax.YData = ["P(b)"];

% Plot state distributions.
figure(2)
% Construct 2D version of environment to visualize probabilities.
grid = zeros(11,7);
grid(5,4) = P_x1(1);
grid(6,5) = P_x1(2);
grid(7,4) = P_x1(3);
grid(6,3) = P_x1(4);
heatmap(grid, 'ColorLimits',[0 1], 'Colormap', cool, 'FontSize', 5);
ax = gca;
ax.Title = "P(x1|x0)";
ax.XData = ["x=-3" "-2" "-1" "0" "1" "2" "3"];
ax.YData = ["y=5" "4" "3" "2" "1" "0" "-1" "-2" "-3" "-4" "-5"];
end

%% Update posterior over beta given measured u0
%  Apply Bayesian update:  
%       P'(beta | x0, u0) \propto P(u0 | x0; beta)*P(beta)
function P_b_new = update_Pb(u0, P_b, x0, g, betas)

denominator = Pu_given_x_b(u0, x0, g, betas(1))*P_b(1) + ...
                    Pu_given_x_b(u0, x0, g, betas(2))*P_b(2) + ...
                    Pu_given_x_b(u0, x0, g, betas(3))*P_b(3) + ...
                    Pu_given_x_b(u0, x0, g, betas(4))*P_b(4);

P_b_new = [0,0,0,0];
for bidx=1:length(betas)
    numerator = Pu_given_x_b(u0, x0, g, betas(bidx))*P_b(bidx);
    P_b_new(bidx) = numerator/denominator;
end

end

%% Update state distribution.
%  Marginalize over betas and controls to get state distribution at next
%  timestep:
%       P(x1 | x0) = \sum_beta (\sum_u0 P(x1 | x0, u0)*P(u0 | x0; beta))*P(beta)
function P_x1 = update_Px(x0, g, betas, P_b)

P_x1 = [0,0,0,0]; 
for u=1:4
    p = 0;
    for bidx=1:length(betas)
        p = p + Pu_given_x_b(u, x0, g, betas(bidx))*P_b(bidx);
    end
    P_x1(u) = p;
end

end

%% Compute the probability of a specific action given a state and beta value.
%  Boltzmann observation model:
%       P(u | x0; beta) \propto e^{-beta*Q(x0, u)}
%  where the Q-function is simply:
%       Q(x0, u) = ||u||_2 + ||x0 + u - g||_2
%  and we assume that each control action is norm 1. 
function prob = Pu_given_x_b(u, x0, g, beta)
    global UP DOWN LEFT RIGHT
    
    % Compute the next states that we could possibly get to given our
    % dynamics. 
    xu1 = dynamics(x0,UP);
    xu2 = dynamics(x0,RIGHT);
    xu3 = dynamics(x0,DOWN);
    xu4 = dynamics(x0,LEFT);
    
    % Compute the Q-value of each possible state-action pair.
    Qx0u1 = 1 + norm(xu1 - g);
    Qx0u2 = 1 + norm(xu2 - g);
    Qx0u3 = 1 + norm(xu3 - g);
    Qx0u4 = 1 + norm(xu4 - g);

    % Compute the denominator by summing over all possible actions.
    normalizer = exp(-beta * Qx0u1) + exp(-beta * Qx0u2) + ...
                    exp(-beta * Qx0u3) + exp(-beta * Qx0u4);

    % Compute the numerator by evaluating the likelihood of the given action. 
    if u == UP 
        prob = exp(-beta * Qx0u1)/normalizer;
    elseif u == RIGHT 
        prob = exp(-beta * Qx0u2)/normalizer;
    elseif u == DOWN 
        prob = exp(-beta * Qx0u3)/normalizer;
    elseif u == LEFT 
        prob = exp(-beta * Qx0u4)/normalizer;
    end
end

%% Dynamics function gives next state we can get to given current state.
function xnext = dynamics(x0,u)
    global UP DOWN LEFT RIGHT
    if u == UP
        xnext = [x0(1); x0(2)+1];
    elseif u == RIGHT
        xnext = [x0(1)+1; x0(2)];
    elseif u == DOWN
        xnext = [x0(1); x0(2)-1];
    elseif u == LEFT
        xnext = [x0(1)-1; x0(2)];
    end
end