%% Clear old variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = frsTwoGoalPrecompute();

%% Create human predictor.
predictor = HumanPredictor(params);

%% Discretize x,y, and probability of beta
xvalues = linspace(params.lowEnv(1), params.upEnv(1), 17);
yvalues = linspace(params.lowEnv(2), params.upEnv(2), 17);
Pbetas = linspace(0,1,11);

%% File structure for saving
frsDir = '/home/ext_drive/somilb/safe_navigation_data/safe_hallucination/frs_two_goal/';

for pb = Pbetas
    for xval = xvalues
        for yval = yvalues
            % Set the new initial condition
            initz = [xval; yval; pb];

            predictor.human.x = initz;
            predictor.human.xhist = initz;
            predictor.zcurr = initz;

            % Compute the predictions
            predictor.updatePredictions();
            [preds, times] = predictor.getPredictions();

            % Project to 2D.
            [grid2D, preds2D] = proj(params.predGrid, preds, [0,0,1], 'min');

            % Convert to binary map.
            occuMap2D = preds2D;
            occuMap2D = 1*(occuMap2D <= 0) + 0*(occuMap2D > 0);

            % Filenames and paths.
            betaFolder = strcat('pb', sprintf('%0.2f', pb));
            predPath = strcat(frsDir, betaFolder);
            xyname = strcat('xval_', sprintf('%0.2f', xval), '_yval_', sprintf('%0.2f', yval));
            occuMat = 'occuMap.mat';
            fullFilePath = strcat(predPath, '/', xyname, '_', occuMat);
            
            % Create directory if it doesnt exist. 
            if ~exist(predPath, 'dir')
              mkdir(predPath)
            end

            % Save full horizon FRS projected into 2D and masked.
            predictions = occuMap2D;
            predGrid = grid2D;
            predTmin = params.tMin;
            predTmax = params.tMax;
            predDt = params.dt;
            predTimes = times;
            save(fullFilePath, 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');

        end
    end
end