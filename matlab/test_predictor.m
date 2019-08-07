%% Setup. 
clear all
clf

% Initial state of the human (in cells).
xinit = [8;5];

% Known goal (in cells).
goal = [2;5];

% Num grid cells in X and Y dimension. 
gridDims = [10,10];

% Set of values beta can take.
betas = [1, 10, 100, 1000];

% Prior over beta.
%P_betas = [0.1 ,0.1, 0.1, 0.7]; 
%P_betas = [1,0,0,0]; % prior that gives us FRS.
P_betas = [0,0,0.1,0.9];

%% Construct predictor.
predictor = Predictor(xinit, goal, gridDims, betas, P_betas);

H = 5; % prediction horizon
ctrlThresh = 0.01;
predGrids = predictor.predictDeterministic(xinit, H, ctrlThresh);

% Plot beta distribution and the state preds.
predictor.plot(predGrids, 1:H+1);