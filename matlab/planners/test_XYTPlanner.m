% Planner discretization parameters.
xDisc = 0.25;
yDisc = 0.25;
tDisc = 0.001;
edgeVel = 0.5;

heurWeight = 5;

planner = XYTPlanner(xDisc, yDisc, tDisc, heurWeight, edgeVel);

% Set up the state bounds.
planner.stateBounds('x') = [-0.1, 3];
planner.stateBounds('y') = [-0.1, 3];
planner.stateBounds('t') = [0, 12]; % Planning horizon.

% Set up the occupancy grid.
xBounds = planner.stateBounds('x');
yBounds = planner.stateBounds('y');
tBounds = planner.stateBounds('t');
planner.staticObsMap = OccupancyGrid(xDisc, yDisc, tDisc, ...
                                     xBounds(1), xBounds(2), ...
                                     yBounds(1), yBounds(2), ...
                                     tBounds(1), tBounds(2));

obsBoundsX = [0.5, 1.5];
obsBoundsY = [0.75, 1.5];
% for x = obsBoundsX(1):xDisc:obsBoundsX(2)
%     for y = obsBoundsY(1):yDisc:obsBoundsY(2)
%         planner.staticObsMap.setData(x, y, 0, 1);
%     end
% end

% dynObsStartXY = [0; 2.5];
dynObsStartXY = [0.3; 1];
% dynObsVel = [0.25; -0.1];
dynObsVel = [0.2; 0];
dynObsDim = [0.75; 0.75];

planner.dynObsMap = OccupancyGrid(xDisc, yDisc, tDisc, ...
                                  xBounds(1), xBounds(2), ...
                                  yBounds(1), yBounds(2), ...
                                  tBounds(1), tBounds(2));
for t = tBounds(1):tDisc:tBounds(2)
    lowerBounds = dynObsStartXY + t * dynObsVel;
    upperBounds = dynObsStartXY + t * dynObsVel + dynObsDim;

    k = planner.dynObsMap.timeToIndex(t);
    planner.dynObsMap.data{k} = zeros(ceil((planner.dynObsMap.xMax - ...
                                            planner.dynObsMap.xMin) / ...
                                           planner.dynObsMap.xDisc), ...
                                      ceil((planner.dynObsMap.yMax - ...
                                            planner.dynObsMap.yMin) / ...
                                           planner.dynObsMap.yDisc));

    for x = lowerBounds(1):xDisc:upperBounds(1)
        for y = lowerBounds(2):yDisc:upperBounds(2)
            planner.dynObsMap.setData(x, y, t, 1);
        end
    end

    % for x = xBounds(1):xDisc:xBounds(2)
    %     for y = yBounds(1):yDisc:yBounds(2)
    %         planner.dynObsMap.setData(x, y, t, 0);
    %     end
    % end
end

goalXY = [1.8; 1.8];
goalTol = 0.2;

% Plan the trajectory.
startStateCont = [0; 0; 0];
traj = planner.plan(startStateCont, goalXY, goalTol);

figure;
hold on;
axis equal;

% Draw the goal and goal tolerance.
scatter([goalXY(1)], [goalXY(2)]);
drawCircle(goalXY, goalTol);

% Draw the obstacle.
% obs = fill([obsBoundsX(1), obsBoundsX(1), obsBoundsX(2), obsBoundsX(2)], ...
%            [obsBoundsY(1), obsBoundsY(2), obsBoundsY(2), obsBoundsY(1)], ...
%            'r');
% set(obs, 'facealpha', 0.5);

lowerBounds = dynObsStartXY;
upperBounds = dynObsStartXY + dynObsDim;

dynObs = fill([lowerBounds(1), lowerBounds(1), upperBounds(1), upperBounds(1)], ...
              [lowerBounds(2), upperBounds(2), upperBounds(2), lowerBounds(2)], ...
              'g');
dynObs.set('facealpha', 0.5);

traj.draw(false);

lastT = traj.contStates{end}(3);
samplesX = [];
samplesY = [];
t = 0;
dt = 0.01;
state = traj.getState(t);
stateVis = drawCircle(state(1:2), 0.05);

while t < lastT
    state = traj.getState(t);

    t = t + dt;
    samplesX = [samplesX, state(1)];
    samplesY = [samplesY, state(2)];

    % Update the state.
    stateVis = drawCircle(state(1:2), 0.05);

    % Update the dynmaic obstacle.
    lowerBounds = dynObsStartXY + t * dynObsVel;
    upperBounds = dynObsStartXY + t * dynObsVel + dynObsDim;
    dynObs.set('XData', [lowerBounds(1), lowerBounds(1), upperBounds(1), ...
                        upperBounds(1)]);
    dynObs.set('YData', [lowerBounds(2), upperBounds(2), upperBounds(2), ...
                        lowerBounds(2)]);

    drawnow;
    pause(dt);
end

scatter(samplesX, samplesY);

function circ = drawCircle(origin, radius)
    th = 0:pi/50:2*pi;
    xunit = radius * cos(th) + origin(1);
    yunit = radius * sin(th) + origin(2);
    circ = plot(xunit, yunit);
end
