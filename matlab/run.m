%% Clear old figure plotting and variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = dubinsCarGaussianHuman();

%% Create human predictor.
predictor = HumanPredictor(params);

%% Create robot motion planner.
planner = LatticePlanner(params.xDisc, params.yDisc, ...
                        params.thetaDisc, params.vDisc, ...
                        params.tDisc, params.heurWeight);
                    
% Set up the state bounds.
planner.stateBounds('x') = [params.lowEnv(1), params.upEnv(1)];
planner.stateBounds('y') = [params.lowEnv(2), params.upEnv(2)];
planner.stateBounds('v') = [0.1, 3];
planner.stateBounds('theta') = [-pi, pi];
planner.stateBounds('t') = [0, 15]; % Planning horizon.
                  
% Create the static obstacle map.
planner.staticObsMap = OccupancyGrid(params.xDisc, params.yDisc, params.tDisc, ...
                                     params.lowEnv(1), params.upEnv(1), ...
                                     params.lowEnv(2), params.upEnv(2), ...
                                     params.tMin, params.tMax);
% Create the dynamic obstacle map.
tDisc = (params.tMax - params.tMin)/params.dt;
planner.dynObsMap = OccupancyGrid(params.N(1), params.N(2), tDisc, ...
                                  params.lowEnv(1), params.upEnv(1), ...
                                  params.lowEnv(2), params.upEnv(2), ...
                                  params.tMin, params.tMax);

%% Simulation loop.
hold on

% Create the first prediction.
predictor.updatePredictions();
[preds, times] = predictor.getPredictions();
planner.dynObsMap.fromValueFuns(params.predGrid, preds, times, 0);

% Create the first plan. 
traj = planner.plan(params.robotx0, params.robotGoalXY, params.goalTol);

for t=1:params.T
    
    % Get most recent measurement of where the person is and what action
    % they applied.
    [xHnext, uHcurr] = params.simHuman.simulateAction(params.simDt);
    
    plot(xHnext(1), xHnext(2),'k.', 'MarkerSize', 10);
    xlim([-1,3]);
    ylim([-1,3]);
    pause(0.1);
    
    % Prediction step. 
    predictor.updateState(xHnext, 0);
    predictor.xcurr
    %predictor.updatePredictions();
    
    % Planning step. 
    % [preds, times] = predictor.getPredictions();
    % planner.updatePredictions(preds, times);
    % plan = planner.replan();
    
    % Move the robot.
    % robot.executePlan(plan);
end