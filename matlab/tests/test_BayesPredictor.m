%% This file tests the FULLY BAYESIAN predictor
clf
clear all

% Variance on Gaussian observation model.
sigma1 = pi/8;
sigma2 = pi/8;

% Known goal locations (in m). 
goals = {[1, 1], [1, -1]}; %{[1, tan(sigma1)], [1, -tan(sigma2)]};

% Grid structure definitions
gridMin = [-2,-2];          % Lower & upper bounds of grid (in m)
gridMax = [2,2];
gridDims = [41,41];         % Num grid cells in X and Y dimension. 

% Set the prior over goal 1 and goal 2.
prior = [0.5, 0.5];

% Create the predictor. 
predictor = BayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, gridDims);
%predictor.Pu_given_x_g(2, [21,21], 2);

% --- debugging ------ %
%predictor.plot_Pu_given_x();
% -------------------- %

% Initial human state (in m).
x0 = [0,0];

% Prediction horizon. 
T = 3;                      % horizon in (seconds)
dt = 0.5;  
H = T/dt;                   % horizon in (timesteps)

% Predict!
preds = predictor.predict(x0, H);

gString = createGrid(gridMin, gridMax, gridDims);

% Plot
figure(2)
set(gcf,'color','w');

epsilon = 0;
for t=1:H+1
    p = preds{t};
    sum(sum(p))
    
    % For visualizingi all distributions
    %pcolor(gString.xs{2}, gString.xs{1}, p);
    %colorbar
    %caxis([0,1])
    
    % For visualizing thresholded. 
    p1 = (p > epsilon);
    [M,c1] = contour(gString.xs{2}, gString.xs{1}, p1, [0,0.1]);
    c1.LineWidth = 2;
    c1.Color = 'b';
    %hold on
%     p2 = (p > 0);
%     [M,c2] = contour(gString.xs{2}, gString.xs{1}, p2, [0,0.1]);
%     c2.LineWidth = 2;
%     c2.Color = 'r';
    grid on
    
    title(strcat('t=', num2str((t-1)*dt)));
end