%% Clear old figure plotting and variables.
clf 
clc 
clear all

%% Load the experimental setup.
%params = scenario1();
% params = scenario2();
params = scenario1_splines();

% Load the predictions.
% load('/home/abajcsy/hybrid_ws/src/hallucinate/matlab/data/fixed_human_ours_p05.mat');
% load('/home/abajcsy/hybrid_ws/src/hallucinate/matlab/data/fixed_human_rss_p05.mat');

load(['/home/eratner/Documents/hallucinate/matlab/data/fixed_human_ours_p05.mat']);
% load(['/home/eratner/Documents/hallucinate/matlab/data/fixed_human_rss_p05.mat']);


%% Create human predictor.
predictor = HumanPredictor(params);

%% Create robot motion planner.
planner = LatticePlanner(params.xDisc, params.yDisc, ...
                         params.thetaDisc, params.vDisc, ...
                         params.tDisc, params.heurWeight, ...
                         params.minEdgeTimeSteps, params.maxEdgeTimeSteps);

% Set up the state bounds.
planner.stateBounds('x') = params.xBounds;
planner.stateBounds('y') = params.yBounds;
planner.stateBounds('v') = params.vBounds;
planner.stateBounds('theta') = params.thetaBounds;
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
% TODO (Ellis): This may not work anymore
for t=params.tMin:params.tDisc:params.tMax
    for r = params.staticObsBounds
        rect = r{1};
        planner.staticObsMap.setRectangularObs(rect{1}, rect{2}, t);
    end
end

% Setup the dynamic obstacle map.
predGrid2D = proj(params.predGrid, params.predGrid.xs, [0,0,1]);
planner.dynObsMap = OccupancyGrid(params.xDisc, params.yDisc, params.dt, ...
                                  params.lowEnv(1), params.upEnv(1), ...
                                  params.lowEnv(2), params.upEnv(2), ...
                                  params.tMin, params.tMax, predGrid2D, ...
                                  params.pathToFRSDir);

%% Simulation loop.
hold on

% Initialize the planned trajectory.
contStates = {};
traj = Trajectory(contStates);

% Initialize the robot state.
xRActual = params.simRobot.x; 
xR = params.xR0;
xRLast = xR;

% Get the initial state of the simulated human.
xHnext = params.xH0(1:2);

hh = [];
rh = [];
arh = [];
th = [];

staticMapHandle = [];
dynMapHandle = [];

% pIdx = 2;
for t=0:params.T-1
    fprintf('At time step %d, robot is at (%f, %f, %f, %f)\n', ...
            t, xR(1), xR(2), xR(3), xR(4));
    
    xRLast = xR;
    % Update the state based on the planned trajectory / controls.
    if ~traj.isEmpty()
        if ~params.trajUseControl
            if traj.splineIndex(t * params.simDt) > 0
                xR = traj.getState(t * params.simDt);
            else
                fprintf('(run) Warning: time %f not in planned trajectory\n', ...
                        t * params.simDt);
            end
        else
            % Apply control to robot, and integrate the dynamics to get the
            % next state.
            fprintf('Applying a control...\n');
            tspan = [0 params.simDt];
            [~, soln] = ode45(@(t_, x_) unicycleDynamics(t_, x_, traj.getControl((t - 1) * params.simDt + t_)), ...
                              tspan, ...
                              xR);
            xR = soln(end, :)';
        end
    end

    fprintf('Avg robot speed is %f m/s\n', ...
            norm(xR(1:2) - xRLast(1:2)) / params.simDt);

    % ----- plotting ------ %
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
    pause(0.1*params.simDt);
    % --------------------- %

    % Get most recent measurement of where the person is and what action
    % they applied.
    [xHnext, uHcurr] = params.simHuman.simulate(xHnext, ...
        (t + 1) * params.simDt, ...
        params.simDt);
    
    % Prediction step.
    predictor.updateState(xHnext, uHcurr);
    betaProb = round(predictor.zcurr(end), 1);
    
    % Load the predictions from the appropriate file.
    planner.dynObsMap.loadFromFile(betaProb, params.predHorizon);
    planner.dynObsMap.setHumanState(xHnext);
    
    % Plot the occupancy maps.
    fprintf('Plotting static obs map...\n');
    if isempty(staticMapHandle)
        staticMapHandle = planner.staticObsMap.draw('k');
    else
        staticMapHandle = planner.staticObsMap.draw('k', staticMapHandle);
    end
    
    fprintf('Plotting dyn obs map...\n');
    if isempty(dynMapHandle)
        dynMapHandle = planner.dynObsMap.draw(params.predColor);
    else
        dynMapHandle = planner.dynObsMap.draw(params.predColor, ...
            dynMapHandle);
    end

    % Planning step.
    xRStart = xR;
    tStart = t * params.simDt;
    if ~traj.isEmpty() && traj.splineIndex(tStart) > 0
        xRStart = traj.getState(tStart);
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