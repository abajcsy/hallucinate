%% Clear old variables.
clf 
clc 
clear all

sigmas = [pi/4, pi/4];

goals = {[2, 2], [2, -2]};      % Known goal locations (in m). 

% latticeMin = [-3.1, -3.1];
% latticeMax = [3.1, 3.1];
latticeMin = [-3.1, -1.6];
latticeMax = [3.1, 1.6];

sideLength = 0.1;
v = 0.6;
dt = sideLength / v;
fprintf('Using a time discretization of %f for fixed velocity %f\n', ...
    dt, v);

T = 6;
horizon = floor(T/dt);                % Prediction horizon in (timesteps).

goalPriors = [0.5, 0.5];             % Set the prior over goal 1 and goal 2.

%% Create the predictor. 
predictor = OptimizedLatticeBayesPredictor(latticeMin, latticeMax, ...
                sideLength, goals, sigmas, goalPriors);

%% Discretize x,y, and probability of beta
% xvalues = linspace(latticeMin(1), latticeMax(1), 17);
% yvalues = linspace(latticeMin(2), latticeMax(2), 17);

xvalues = linspace(-3, 3, 13);
yvalues = linspace(-1.5, 1.5, 7);

Pbetas = linspace(0, 1, 11);

%% File structure for saving
frsDir = '/home/eratner/research/icra2020/data/bayes_two_goal/';
% frsDir = '/home/ext_drive/somilb/safe_navigation_data/safe_hallucination/bayes_two_goal/';

for pb = Pbetas
    for xval = xvalues
        for yval = yvalues
            predictor.goalPriors = [pb, 1 - pb];
            
            % Compute the predictions.
            fprintf('Computing predictions for prior %f, x %f, y %f...\n', ...
                    pb, xval, yval);
            tstart = cputime;
            preds = predictor.predict(xval, yval, horizon);
            tend = cputime;
            fprintf('done. Took %f s of CPU time.\n', tend - tstart);

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
            [X, Y] = predictor.getLatticeMeshgrid();
            predGridX = X;
            predGridY = Y;
            predTmin = 0;
            predTmax = T;
            predDt = dt;
            save(fullFilePath, 'predictions', 'predGridX', 'predGridY', ...
                'predDt', 'predTmin', 'predTmax');

        end
    end
end

