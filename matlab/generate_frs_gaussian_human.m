%% Clear old variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = frsPrecompute();

%% Create human predictor.
predictor = HumanPredictor(params);

%% Discretize probability of beta
Pbetas = linspace(0,1,11);

%% File structure for saving
%frsDir = '/home/abajcsy/hybrid_ws/src/hallucinate/matlab/frs/';
frsDir = '/home/andreabajcsy/hybrid_ws/src/hallucinate/matlab/frs/';

for pb = Pbetas
    % Set the new initial condition
    initz = [params.xH0; pb];
   
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
    betaFolder = strcat('pb', num2str(pb));
    predPath = strcat(frsDir, betaFolder);
    
    horiz6 = strcat(predPath, '/h6/');
    horiz4 = strcat(predPath, '/h4/');
    horiz2 = strcat(predPath, '/h2/');
    
    frsMat = 'fullFRS.mat';
    occuMat = 'occuMap.mat';
    
    % Save full horizon FRS.
    predictions = preds;
    predGrid = params.predGrid;
    predTmin = params.tMin;
    predTmax = params.tMax;
    predDt = params.dt;
    predTimes = times;
    
    % 6 second full pred.
    save(strcat(horiz6, frsMat), 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');
    
    % 6 second occupancy map.
    predictions = occuMap2D;
    predGrid = grid2D;
    save(strcat(horiz6, occuMat), 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');
    
    % 4 second full pred.
    idx4sec = 4/params.dt;
    predictions = preds(:,:,:,1:idx4sec);
	predTmax = 4;
    predTimes = times(1:idx4sec);
    predGrid = params.predGrid;
    save(strcat(horiz4, frsMat), 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');
    
    % 4 second occupancy map.
    predictions = occuMap2D(:,:,1:idx4sec);
    predGrid = grid2D;
    save(strcat(horiz4, occuMat), 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');
    
    % 2 second full pred.
    idx2sec = 2/params.dt;
    predictions = preds(:,:,:,1:idx2sec);
    predTmax = 2;
    predTimes = times(1:idx2sec);
    predGrid = params.predGrid;
    save(strcat(horiz2, frsMat), 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');
    
    % 2 second occupancy map.
    predictions = occuMap2D(:,:,1:idx2sec);
    predGrid = grid2D;
    save(strcat(horiz2, occuMat), 'predictions', 'predGrid', 'predTimes', 'predDt', 'predTmin', 'predTmax');
end