%% Clear old variables.
clf 
clc 
clear all

sigma1 = pi/4;                  % Variance on Gaussian observation model.
sigma2 = pi/4;

goals = {[2, 2], [2, -2]};      % Known goal locations (in m). 

gridMin = [-4,-4];              % Lower & upper bounds of grid (in m).
gridMax = [4,4];

r = 0.1;                        % Grid cell size.
v = 0.6;                        % Velocity of human. 
dt = r / v;                     % Prediction timestep. 

T = 2;                          % Prediction horizon in (seconds).
H = floor(T/dt);                % Prediction horizon in (timesteps).

prior = [0.5, 0.5];             % Set the prior over goal 1 and goal 2.

%% Create the predictor. 
predictor = LatticeBayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, r);

%% Discretize x,y, and probability of beta
xvalues = linspace(params.lowEnv(1), params.upEnv(1), 17);
yvalues = linspace(params.lowEnv(2), params.upEnv(2), 17);
Pbetas = linspace(0,1,11);

%% File structure for saving
frsDir = '/home/ext_drive/somilb/safe_navigation_data/safe_hallucination/bayes_two_goal/';

for pb = Pbetas
    for xval = xvalues
        for yval = yvalues
            % Set the new initial state (in m) and prior for predictor. 
            x0 = [xval, yval];
            predictor.prior = containers.Map([1:length(prior)], [pb, 1 - pb]);
            
            % Compute the predictions.
            preds = predictor.predict(x0, H);

            % Filenames and paths.
            betaFolder = strcat('pb', sprintf('%0.2f', pb));
            predPath = strcat(frsDir, betaFolder);
            xyname = strcat('xval_', sprintf('%0.2f', xval), '_yval_', sprintf('%0.2f', yval));
            predMat = 'predMap.mat';
            fullFilePath = strcat(predPath, '/', xyname, '_', predMat);
            
            % Create directory if it doesnt exist. 
            if ~exist(predPath, 'dir')
              mkdir(predPath)
            end

            % Save full horizon FRS projected into 2D and masked.
            predictions = preds;
            predGrid = grid2D;
            predTmin = 0;
            predTmax = T;
            predDt = dt;
            save(fullFilePath, 'predictions', 'predGrid', 'predDt', 'predTmin', 'predTmax');

        end
    end
end

