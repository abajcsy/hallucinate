%% Clear old variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = frsPrecompute();

%% Create human predictor.
predictor = HumanPredictor(params);

%% Discretize probability of beta
Pbetas = linspace(0,1,10);

%% File structure for saving
frsDir = '/home/abajcsy/hybrid_ws/src/hallucinate/matlab/frs/';

for pb = Pbetas
    % Set the new initial condition
    initz = [params.xH0; pb];
    params.humanModel.zcurr = params.z0;
    
    % Compute the predictions
    predictor.updatePredictions();
	[preds, times] = predictor.getPredictions();
    
    % Project to 2D.
    [grid2D, preds2D] = proj(params.predGrid, preds, [0,0,1]);
    
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
    save(strcat(horiz6, frsMat), preds, times);
    save(strcat(horiz6, occuMat), occyMap2D, times);
    
    idx4sec = 4/params.dt;
    save(strcat(horiz4, frsMat), preds(:,:,:,1:idx4sec), times(1:idx4sec));
    save(strcat(horiz4, frsMat), occuMap2D(:,:,1:idx4sec), times(1:idx4sec));
    
    idx2sec = 2/params.dt;
    save(strcat(horiz2, frsMat), preds(:,:,:,1:idx2sec), times(1:idx2sec));
    save(strcat(horiz2, frsMat), occuMap2D(:,:,1:idx2sec), times(1:idx2sec));
end