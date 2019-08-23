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
planner.stateBounds('x') = params.xBounds;
planner.stateBounds('y') = params.yBounds;
planner.stateBounds('v') = params.vBounds;
planner.stateBounds('theta') = params.thetaBounds;
planner.stateBounds('t') = params.timeBounds;
                  
% Create the static obstacle map.
Nx = ceil((params.upEnv(1) - params.lowEnv(1))/params.xDisc);
Ny = ceil((params.upEnv(2) - params.lowEnv(2))/params.yDisc);
staticGrid2D = createGrid(params.lowEnv, params.upEnv, [Nx; Ny]);
planner.staticObsMap = OccupancyGrid(params.xDisc, params.yDisc, params.tDisc, ...
                                     params.lowEnv(1), params.upEnv(1), ...
                                     params.lowEnv(2), params.upEnv(2), ...
                                     params.tMin, params.tMax);
% Setup the static obstacles.
for t=params.tMin:params.tDisc:params.tMax
    for r = params.staticObsBounds
        rect = r{1};
        planner.staticObsMap.setRectangularObs(rect{1}, rect{2}, t);
    end
end
                                 
% Setup the dynamic obstacle map.
tDisc = (params.tMax - params.tMin)/params.dt;
planner.dynObsMap = OccupancyGrid(params.N(1), params.N(2), tDisc, ...
                                  params.lowEnv(1), params.upEnv(1), ...
                                  params.lowEnv(2), params.upEnv(2), ...
                                  params.tMin, params.tMax);

%% Simulation loop.
hold on

% Create the first prediction.
%predictor.updatePredictions();
%[preds, times] = predictor.getPredictions();
%planner.dynObsMap.fromValueFuns(params.predGrid, preds, times, 0);

% Create the first plan. 
traj = planner.plan(params.xR0, params.goalRXY, params.goalTol);

% Apply control to robot
uR = traj.getControl(0);
params.simRobot.updateState(uR, params.simDt, params.simRobot.x);
xR = params.simRobot.x;

% Get the initial state of the simulated human.
xHnext = params.xH0(1:2);

for t=1:params.T
        
    % ----- plotting ------ %
    planner.staticObsMap.draw(staticGrid2D);    % draw the static obstacle
    %planner.dynObsMap.draw(grid2D);      % draw the dynamic obstacle
    plotAgent(xHnext, 'b');         % plot human
    plotAgent(xR, 'k');             % plot robot
    plot(params.goalRXY(1), params.goalRXY(2), 'ko','MarkerSize', 8); % plot robot goal
    traj.draw();                    % plot trajectory
    xlim([params.lowEnv(1),params.upEnv(1)]);
    ylim([params.lowEnv(2),params.upEnv(2)]);
    pause(0.1);
    % --------------------- %
    
    % Get most recent measurement of where the person is and what action
    % they applied.
    [xHnext, uHcurr] = params.simHuman.simulateAction(params.simDt);

    % Prediction step. 
    %predictor.updateState(xHnext, uHcurr);
    %predictor.updatePredictions();
    
    % Planning step. 
    %[preds, times] = predictor.getPredictions();
    %planner.dynObsMap.fromValueFuns(params.predGrid, preds, times, t*params.simDt);
    traj = planner.plan(params.xR0, params.goalRXY, params.goalTol);

    % Apply control to robot.
    uR = traj.getControl(t*params.simDt);
    params.simRobot.updateState(uR, params.simDt, xR);
    xR = params.simRobot.x;

end

%% Plots human or robot.
% Inputs:
%   x [vector]  - 3D/4D state of agent
% Ouput:
%   c   - handle for figure
function c = plotAgent(x, color)
    c = {};
    c{1} = plot(x(1), x(2), 'ko','MarkerSize', 8, ...
        'MarkerEdgeColor', color, 'MarkerFaceColor', color);

%     % Plot heading.
%     center = x(1:2);
% 
%     if length(x) >= 3
%         % Rotation matrix.
%         R = [cos(x(3)) -sin(x(3)); 
%              sin(x(3)) cos(x(3))];
%         % Heading pt.
%         hpt = [0.2; 0];
%         hptRot = R*hpt + center;
%         p2 = plot([center(1) hptRot(1)], [center(2) hptRot(2)], color, 'LineWidth', 1.5);
%         p2.Color(4) = 1.0;
%         c{2} = p2;
%     end
end