%% Clear old figure plotting and variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = scenario1();
%params = scenario2();

% Load the predictions.
load('/home/abajcsy/hybrid_ws/src/hallucinate/matlab/data/fixed_human_ours_p05.mat');
% load('/home/abajcsy/hybrid_ws/src/hallucinate/matlab/data/fixed_human_rss_p05.mat');

% load(['/home/eratner/Documents/hallucinate/matlab/data/' ...
%       'fixed_human_ours_p05.mat']);

%% Create human predictor.
predictor = HumanPredictor(params);

%% Create robot motion planner.
% planner = LatticePlanner(params.xDisc, params.yDisc, ...
%                          params.thetaDisc, params.vDisc, ...
%                          params.tDisc, params.heurWeight, ...
%                          params.minEdgeTimeSteps, params.maxEdgeTimeSteps);

planner = XYTPlanner(params.xDisc, params.yDisc, params.tDisc, ...
                     params.heurWeight, params.deltaTCont);


% Set up the state bounds.
planner.stateBounds('x') = params.xBounds;
planner.stateBounds('y') = params.yBounds;
% planner.stateBounds('v') = params.vBounds;
% planner.stateBounds('theta') = params.thetaBounds;
planner.stateBounds('t') = params.timeBounds;

% Create the static obstacle map.
Nx = round((params.upEnv(1) - params.lowEnv(1))/params.xDisc) + 1;
Ny = round((params.upEnv(2) - params.lowEnv(2))/params.yDisc) + 1;
staticGrid2D = createGrid(params.lowEnv, params.upEnv, [Nx; Ny]);
planner.staticObsMap = OccupancyGrid(params.xDisc, params.yDisc, params.tDisc, ...
                                     params.lowEnv(1), params.upEnv(1), ...
                                     params.lowEnv(2), params.upEnv(2), ...
                                     params.tMin, params.tMax, staticGrid2D);
% Setup the static obstacles.
for t=params.tMin:params.tDisc:params.tMax
    for r = params.staticObsBounds
        rect = r{1};
        planner.staticObsMap.setRectangularObs(rect{1}, rect{2}, t);
    end
end

% Setup the dynamic obstacle map.
tDisc = (params.tMax - params.tMin)/params.dt;
predGrid2D = proj(params.predGrid, params.predGrid.xs, [0,0,1]);
planner.dynObsMap = OccupancyGrid(params.xDisc, params.yDisc, params.tDisc, ...
                                  params.lowEnv(1), params.upEnv(1), ...
                                  params.lowEnv(2), params.upEnv(2), ...
                                  params.tMin, params.tMax, predGrid2D);


%% Simulation loop.
hold on

% predictor.updatePredictions();
% [preds, times] = predictor.getPredictions();
planner.dynObsMap.fromValueFuns(predGrid2D, humanPreds{1}, predsTimes{1}, 0);

% Initialize the planned trajectory.
contStates = {};
% traj = Trajectory(contStates);
traj = XYTTrajectory(contStates);

% Initialize the robot state.
xRActual = params.xR0(1:4, :);
xR = params.xR0(1:2);
xRLast = xR;

% Get the initial state of the simulated human.
xHnext = params.xH0(1:2);

hh = [];
rh = [];
arh = [];
th = [];
pIdx = 2;
for t=0:params.T-1
    xRLast = xR;
    % Update the state based on the planned trajectory / controls.
    if ~traj.isEmpty()
        if ~params.trajUseControl
            if traj.inBounds(t * params.simDt)
                xR = traj.getState(t * params.simDt);
                xR = xR(1:2);
            else
                fprintf('(run) Warning: time %f not in planned trajectory\n', ...
                        t * params.simDt);
            end
        else
            % Apply control to robot, and integrate the dynamics to get the
            % next state.
            fprintf('Applying a control...\n');
            tspan = [0 params.simDt];
            [~, soln] = ode45(@(t_, x_) unicycleDynamics(t_, x_, ...
                getControl((t - 1) * params.simDt + t_, x_, traj)), tspan, xRActual);
            xRActual = soln(end, :)';
            xR = xRActual(1:2);
        end
    end

    fprintf('Robot speed is %f m/s\n', ...
            norm(xR - xRLast) / params.simDt);

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
    end

    if ~isempty(arh)
        delete(arh{1});
        delete(arh{2});
    end

    if ~isempty(th)
        for i=1:length(th)
            delete(th{i});
        end
    end
    th = traj.draw();                   % plot trajectory
    rh = plotAgent(xR, 'r');            % plot robot
    arh = plotAgent(xRActual, 'g');

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
    if ~traj.isEmpty() && traj.inBounds(tStart)
        xRStart = traj.getState(tStart);
        xRStart = xRStart(1:2);
    end

    if traj.isEmpty() || mod(t, params.replanAfterSteps) == 0
        traj = planner.plan([xRStart; tStart], params.goalRXY, ...
                            params.goalTol);
    end
end

%% Gets the 4D Unicycle dynamics.
function dxdt = unicycleDynamics(t, x, u)
    dxdt = [x(4) * cos(x(3));
            x(4) * sin(x(3));
            u(1);
            u(2)];
end

%% Gets the 3D Dubins car dynamics.
function dxdt = dubinsCarDynamics(t, x, u, v)
    dxdt = [v * cos(x(3));
            v * sin(x(3));
            u];
end

%% Simple controller.
function u = getControl(t, x, trajRef)
    kPropLinAcc = 2;
    kPropAngAcc = 0.5;

    xRef = trajRef.getState(t);
    posErr = [xRef(1) - x(1);
              xRef(2) - x(2)];

    % Acceleration is proportional to position error along current
    % direction of travel.
    dir = [cos(x(3)); sin(x(3))];
    a = kPropLinAcc * (posErr' * dir);

    % Angular velocity is proportional to the difference between the
    % current heading and the reference heading.
    headingErr = xRef(3) - x(3);
    omega = kPropAngAcc * headingErr;
    u = [omega; a];
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