%% Testing the optimized version of Bayesian prediction on a equilateral 
%% triangular lattice.
clf
clear all

% Variance on Gaussian observation model.
sigmas = [pi/4, pi/4];

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]};

% Grid structure definitions
latticeMin = [-4, -4];          
latticeMax = [4, 4];

% Set the prior over goal 1 and goal 2.
% goalPriors = [0.5, 0.5];
goalPriors = [0.8, 0.2];

% Equilateral triangle side length.
sideLength = 0.1;

time_s = 3;
v = 0.6;
dt = sideLength/v;

predictionHorizon = ceil(time_s/dt);
x0 = 0;
y0 = 0;

% Create the predictor. 
predictor = OptimizedLatticeBayesPredictor(latticeMin, latticeMax, sideLength, ...
    goals, sigmas, goalPriors);

preds = predictor.predict(x0, y0, predictionHorizon);

% Save the predictions
filename = strcat('./data_for_paper/example1/', 'predictions_bayesian_prior_', num2str(goalPriors(1)), '.mat');
save(filename, 'preds', 'predictor', 'v', 'time_s', 'dt', 'goals', '-v7.3');

