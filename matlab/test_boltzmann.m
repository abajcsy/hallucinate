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

% Environment and gridding information.
realLow = [-3,-5];
realUp = [3,5];
resX = 1;
resY = 1;

% Initial (x,y) state of the human (in meters).
x0 = [0;0];

% Known goal (in meters).
g = [0;4];

% Set of values beta can take.
betas = [1, 10, 100, 1000];
str_betas = {};
for i=1:length(betas)
    str_betas{i} = num2str(betas(i));
end

% Prior over beta.
P_b = [0.1 ,0.1, 0.1, 0.7];

%% Update prior over \beta given new measurement.
% Simulate a bunch of "good" actions (UP) and "bad" actions (DOWN).
P_b_new = P_b;

num_good_us = 10;
num_bad_us = 0;

for i=1:num_good_us
    P_b_new = update_Pb(UP, P_b_new, x0, g, betas);
end

for i=1:num_bad_us
    P_b_new = update_Pb(RIGHT, P_b_new, x0, g, betas);
end

%% Predict next state distribution given your posterior.
H = 5; % prediction horizon
P_x1 = predict(x0, g, H, P_b_new, betas);

%% Plotting.
plot(P_x1, P_b_new, str_betas)
end

%% Predicts the state distribution H steps into the future given x0 
%  Computes:
%       P(xt | x0) = \sum_xt-1 P(x_t | x_t-1)*P(x_t-1 | x_0)
%                  
%  Given:
%       x0    -- initial state
%       g     -- goal location
%       H     -- prediction horizon
%       P_b   -- prior over beta
%       betas -- values that beta can take
%  Output:
%       preds -- cell array indexed by 1:H with corresponding state
%                distributions
function preds = predict(x0, g, H, P_b, betas)

%     assume P(xt | x0) = 0 for all xt
%     for t in {1,2,...H}
%        for xt in X
%            for beta in B
%                for ut-1 in U
%                    xt-1 = finv(xt, ut-1)
%                    p(xt | x0) += P(beta)*P(ut-1|xt-1,beta)*P(xt-1|x0)


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

    % Normalization trick to improve numerical stability.
    % See: http://cs231n.github.io/linear-classify/#softmax
    offset = 0; %max([Qx0u1, Qx0u2, Qx0u3, Qx0u4]);
    
    % Compute the denominator by summing over all possible actions.
    normalizer = exp(-log(beta) * Qx0u1 - offset) + exp(-log(beta) * Qx0u2 - offset) + ...
                    exp(-log(beta) * Qx0u3 - offset) + exp(-log(beta) * Qx0u4 - offset);

    % Compute the numerator by evaluating the likelihood of the given action. 
    if u == UP 
        prob = exp(-log(beta) * Qx0u1 - offset)/normalizer;
    elseif u == RIGHT 
        prob = exp(-log(beta) * Qx0u2 - offset)/normalizer;
    elseif u == DOWN 
        prob = exp(-log(beta) * Qx0u3 - offset)/normalizer;
    elseif u == LEFT 
        prob = exp(-log(beta) * Qx0u4 - offset)/normalizer;
    end
end

%% Gets the 1-step reachable states from x0
function x1 = oneStepReachableStates(x0)
    global UP DOWN LEFT RIGHT
    xu1 = dynamics(x0,UP);
    xu2 = dynamics(x0,RIGHT);
    xu3 = dynamics(x0,DOWN);
    xu4 = dynamics(x0,LEFT);
    x1 = [xu1, xu2, xu3, xu4];
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

%% Plotting. 
function plot(P_x1, P_b, str_betas)
    % Plot the posterior over beta.
    figure(1)
    heatmap(P_b, 'ColorLimits',[0 1]);
    ax = gca;
    ax.XData = str_betas; 
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