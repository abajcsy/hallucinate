%% Clear old figure plotting and variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = dubinsCarGaussianHuman();

%% Create human predictor.
predictor = HumanPredictor(params);

%% Create robot motion planner.
%   TBD. 

%% Simulation loop.
hold on
    
% Store the initial prior over beta = 0.
pbeta = params.x0(3);
for t=1:params.T
    
    % Get most recent measurement of where the person is and what action
    % they applied.
    [xHnext, uHcurr] = params.simHuman.simulateAction(params.simDt);
    
    plot(xHnext(1), xHnext(2),'k.', 'MarkerSize', 10);
    xlim([-1,3]);
    ylim([-1,3]);
    pause(0.01);
    
    % Prediction step. 
    predictor.updateState([xHnext; pbeta], uHcurr);
    %predictor.updatePredictions();
    
    % Planning step. 
    % [preds, times] = predictor.getPredictions();
    % planner.updatePredictions(preds, times);
    % plan = planner.replan();
    
    % Move the robot.
    % robot.executePlan(plan);
end