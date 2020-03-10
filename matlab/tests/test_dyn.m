clf
clear all

% Variance on Gaussian observation model.
sigma1 = pi/4;
sigma2 = pi/4;

% Known goal locations (in m). 
goals = {[4, 4], [4, 0]};

% Grid structure definitions
gridMin = [0,0];          % Lower & upper bounds of grid (in m)
gridMax = [4,4];

% Set the prior over goal 1 and goal 2.
prior = [0.9, 0.1];

% Grid cell size.
r = 1;

% HMM model parameter.
hmmParam = 0.6;

% Create the predictor. 
predictor = DynamicLatticeBayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, r, hmmParam);

% Predict!
x0 = [0,2];
H = 3;
preds = predictor.predict(x0, H);