% xDisc = 0.25;
% yDisc = 0.25;
xDisc = 0.15;
yDisc = 0.15;
thetaDisc = pi / 45;
vDisc = 0.5;
tDisc = 0.25;

% heurWeight = 1.25;
heurWeight = 3;

planner = LatticePlanner(xDisc, yDisc, thetaDisc, vDisc, tDisc, heurWeight);

% Set up the state bounds.
planner.stateBounds('x') = [-0.1, 3];
planner.stateBounds('y') = [-0.1, 3];
planner.stateBounds('v') = [0.1, 3];
planner.stateBounds('theta') = [-pi, pi];
planner.stateBounds('t') = [0, 15]; % Planning horizon.

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
dynObsVel = [0.75; 0];
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

contVal = 1.25;
fprintf("cont = %f, disc = %f\n", contVal, planner.contToDisc(contVal, 5));

startState = LatticePlannerState(0, 0, 0, 1, 0);
% startState = LatticePlannerState(2, 2, 22, 1, 0);
[succs, edges, costs] = planner.expand(startState);
fprintf("expanded %d successors\n", length(succs));

% figure;
% hold on;

% for idx = 1:length(succs)
%   % figure;
%   xfunc = @(t) edges{idx}(1, 1) .* t.^3 + edges{idx}(1, 2) .* t.^2 + ...
%           edges{idx}(1, 3) .* t + edges{idx}(1, 4);
%   yfunc = @(t) edges{idx}(2, 1) .* t.^3 + edges{idx}(2, 2) .* t.^2 + ...
%           edges{idx}(2, 3) .* t + edges{idx}(2, 4);

%   fplot(xfunc, yfunc, [0 tDisc]);
% end

% hold off;

% goalXY = [3; 0];
% goalXY = [1; 1];
% goalXY = [0.5; 2.5];
% goalXY = [1.5; 2];
% goalXY = [0.2; 2.8];
% goalXY = [1.8; 1.25];
goalXY = [1.8; 1.8];
goalTol = 0.2;

% Plan the trajectory.
% startStateCont = [0; 0; 0; 0.5; 0];
startStateCont = [0; 0; 0; 1; 0];
% startStateCont = [0; 0; pi / 8; 1; 0];

traj = planner.plan(startStateCont, goalXY, goalTol);

figure;
hold on;
axis equal;

drawTriangle([startStateCont(1); startStateCont(2)], startStateCont(3), 0.1);
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
              'g')
dynObs.set('facealpha', 0.5);

traj.draw(false);

lastT = traj.contStates{end}(5);
samplesX = [];
samplesY = [];
t = 0;
dt = 0.05;
s = traj.getState(t);
stateVis = drawTriangle([s(1); s(2)], s(3), 0.05);

while t < lastT
    s = traj.getState(t);

    t = t + dt;
    samplesX = [samplesX, s(1)];
    samplesY = [samplesY, s(2)];

    drawTriangle([s(1); s(2)], s(3), 0.05, stateVis);

    % Update the dynmaic obstacle.
    lowerBounds = dynObsStartXY + t * dynObsVel;
    upperBounds = dynObsStartXY + t * dynObsVel + dynObsDim;
    dynObs.set('XData', [lowerBounds(1), lowerBounds(1), upperBounds(1), ...
                        upperBounds(1)]);
    dynObs.set('YData', [lowerBounds(2), upperBounds(2), upperBounds(2), ...
                        lowerBounds(2)]);

    drawnow;
    pause(dt);

    % fprintf("control at %f: \n", t);
    % traj.getControl(t)
end

scatter(samplesX, samplesY);

if length(traj.splines) > 0
    s = traj.contStates{end};
    drawTriangle([s(1); s(2)], s(3), 0.1);
end

function tri = drawTriangle(origin, rot, sideLength, lastTri)
    aLocal = [2*sideLength; 0; 1];
    bLocal = [0; -sideLength / sqrt(2); 1];
    cLocal = [0; sideLength / sqrt(2); 1];

    T = [cos(rot), -sin(rot), origin(1);
         sin(rot), cos(rot), origin(2);
         0, 0, 1];
    a = T * aLocal;
    b = T * bLocal;
    c = T * cLocal;

    if nargin > 3
        lastTri.set('XData', [a(1), b(1), c(1)]);
        lastTri.set('YData', [a(2), b(2), c(2)]);
        tri = lastTri;
    else
        tri = fill([a(1), b(1), c(1)], ...
                   [a(2), b(2), c(2)], ...
                   'b');
        set(tri, 'facealpha', 0.5);
    end
end

function circ = drawCircle(origin, radius)
    th = 0:pi/50:2*pi;
    xunit = radius * cos(th) + origin(1);
    yunit = radius * sin(th) + origin(2);
    circ = plot(xunit, yunit);
end
