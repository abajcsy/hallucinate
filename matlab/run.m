%% Clear old figure plotting and variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = scenario1();
%params = scenario2();

% Load the predictions.
load('/home/abajcsy/hybrid_ws/src/hallucinate/matlab/data/fixed_human_ours_p05.mat');

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
predGrid2D = proj(params.predGrid, params.predGrid.xs, [0,0,1]);
                              
%% Simulation loop.
hold on

% predictor.updatePredictions();
% [preds, times] = predictor.getPredictions();
planner.dynObsMap.fromValueFuns(predGrid2D, humanPreds{1}, predsTimes{1}, 0);

% Initialize the planned trajectory.
contStates = {};
traj = Trajectory(contStates);

% Initialize the robot state.
xR = params.xR0(1:4, :);

% Get the initial state of the simulated human.
xHnext = params.xH0(1:2);

hh = [];
rh = [];
th = [];
pIdx = 2;
for t=0:params.T-1
    % Update the state based on the planned trajectory / controls.
    if ~traj.isEmpty()
        if ~params.trajUseControl
            xR = traj.getState(t * params.simDt);
        else
            % Apply control to robot, and integrate the dynamics to get the next state.
            tspan = [0 params.simDt];
            [~, soln] = ode45(@(t_, x_) unicycleDynamics(t_, x_, ...
                traj.getControl((t - 1) * params.simDt + t_)), tspan, xR);
            xR = soln(end, :)';
        end
    end

    % ----- plotting ------ %
    planner.staticObsMap.draw(staticGrid2D, 'k');           % draw the static obstacle
    planner.dynObsMap.draw(predGrid2D, params.predColor);   % draw the dynamic obstacle
	
    % plot the robot goal and goal tolerance.
    info = [params.goalRXY(1)-params.goalTol params.goalRXY(2)-params.goalTol params.goalTol*2 params.goalTol*2];
    rectangle('Position',info,'Curvature',1, ...
        'FaceColor', [255, 222, 222]/255., 'EdgeColor', [255, 222, 222]/255.);
   	plot(params.goalRXY(1), params.goalRXY(2), ... 
            'ro','MarkerSize', 8);          % plot robot goal
    
    if ~isempty(hh)
        delete(hh{1});
    end
    hh = plotAgent(xHnext, 'b');        % plot human
    
    if ~isempty(rh)
        delete(rh{1});
        delete(rh{2});
    end
    
    if ~isempty(th)
        for i=1:length(th)
            delete(th{i});
        end
    end
    th = traj.draw();                   % plot trajectory
    rh = plotAgent(xR, 'r');            % plot robot
    
    xlim([params.lowEnv(1),params.upEnv(1)]);
    ylim([params.lowEnv(2),params.upEnv(2)]);
    pause(0.1);
    % --------------------- %

    % Get most recent measurement of where the person is and what action
    % they applied.
    [xHnext, uHcurr] = params.simHuman.simulate(xHnext, t*params.simDt, params.simDt);

    % Prediction step.
    %predictor.updateState(xHnext, uHcurr);
    %predictor.updatePredictions();

    % Planning step.
    %[preds, times] = predictor.getPredictions();
    planner.dynObsMap.fromValueFuns(predGrid2D, humanPreds{pIdx}, ...
                predsTimes{pIdx}, ...
                t*params.simDt);
	pIdx = pIdx + 1;

    xRStart = xR;
    tStart = t * params.simDt;
    if ~traj.isEmpty()
        xRStart = traj.getState(tStart);
    end

    if traj.isEmpty() || mod(t, params.replanAfterSteps) == 0
        traj = planner.plan([xRStart; tStart], params.goalRXY, ...
                            params.goalTol);
    end
end

%% Gets the 4D Unicycle dynamics.
function dxdt = unicycleDynamics(t, x, u)
    x(4) = min(x(4), 0.6)
    dxdt = [x(4) * cos(x(3));
            x(4) * sin(x(3));
            u(1);
            u(2)];
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

    % Plot heading.
    center = x(1:2);

    if length(x) >= 3
        % Rotation matrix.
        R = [cos(x(3)) -sin(x(3)); 
             sin(x(3)) cos(x(3))];
        % Heading pt.
        hpt = [0.2; 0];
        hptRot = R*hpt + center;
        p2 = plot([center(1) hptRot(1)], [center(2) hptRot(2)], color, 'LineWidth', 1.5);
        p2.Color(4) = 1.0;
        c{2} = p2;
    end
end