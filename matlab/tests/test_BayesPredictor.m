%% This file tests the FULLY BAYESIAN predictor

% Variance on Gaussian observation model.
sigma1 = pi/8;
sigma2 = pi/8;

% Known goal locations (in m). 
goals = {[1, 1], [1, -1]}; %{[1, tan(sigma1)], [1, -tan(sigma2)]};

% Grid structure definitions
gridMin = [-2,-2];          % Lower & upper bounds of grid (in m)
gridMax = [2,2];
gridDims = [21,21];         % Num grid cells in X and Y dimension. 

% Set the prior over goal 1 and goal 2.
prior = [0.5, 0.5];

% Create the predictor. 
predictor = BayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, gridDims);

% Initial human state (in m).
x0 = [0,0];

% Prediction horizon. 
T = 3;                      % horizon in (seconds)
dt = 0.1;  
H = T/dt;                   % horizon in (timesteps)

% Predict!
preds = predictor.predict(x0, H);

gString = createGrid(gridMin, gridMax, gridDims);

% Plot
figure(2)
for t=1:H
    p = preds{t};
    sum(sum(p))
    pcolor(gString.xs{2}, gString.xs{1}, p);
    colorbar
    caxis([0,1])
end