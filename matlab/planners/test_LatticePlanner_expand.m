% xDisc = 0.1;
% yDisc = 0.1;
xDisc = 0.15;
yDisc = 0.15;
thetaDisc = pi / 45;
vDisc = 0.5;
tDisc = 0.5;

heurWeight = 1.25;
minEdgeTimeSteps = 2;
maxEdgeTimeSteps = 2;

planner = LatticePlanner(xDisc, yDisc, thetaDisc, vDisc, tDisc, heurWeight, ...
                         minEdgeTimeSteps, maxEdgeTimeSteps);

% Set up the state bounds.
planner.stateBounds('x') = [-0.1, 3];
planner.stateBounds('y') = [-0.1, 3];
planner.stateBounds('v') = [0, 0.6];
planner.stateBounds('theta') = [-pi, pi];
planner.stateBounds('t') = [0, 15]; % Planning horizon.

% Set up the (static / dynamic) occupancy grids.
xBounds = planner.stateBounds('x');
yBounds = planner.stateBounds('y');
tBounds = planner.stateBounds('t');
% planner.staticObsMap = OccupancyGrid(xDisc, yDisc, tDisc, ...
%                                      xBounds(1), xBounds(2), ...
%                                      yBounds(1), yBounds(2), ...
%                                      tBounds(1), tBounds(2));
% planner.dynObsMap = OccupancyGrid(xDisc, yDisc, tDisc, ...
%                                   xBounds(1), xBounds(2), ...
%                                   yBounds(1), yBounds(2), ...
%                                   tBounds(1), tBounds(2));

startState = LatticePlannerState(1, 1, 0, 1, 0);
[succs, edges, costs] = planner.expand(startState);
fprintf("expanded %d successors\n", length(succs));

figure;
hold on;

% Draw the start state.
drawTriangle([planner.discToCont(startState.x, 1), ...
              planner.discToCont(startState.y, 2)], ...
             planner.discToCont(startState.theta, 3), 0.01);
         
for idx = 1:length(succs)
    % Draw the successor.
    succState = succs{idx};
    drawTriangle([planner.discToCont(succState.x, 1), ...
                  planner.discToCont(succState.y, 2)], ...
                 planner.discToCont(succState.theta, 3), 0.01);
             
    T = planner.discToCont(succState.t - startState.t, 5);
         
    % Draw the edge.
    a = [edges{idx}(1, 1); edges{idx}(2, 1); edges{idx}(3, 1)];
    b = [edges{idx}(1, 2); edges{idx}(2, 2); edges{idx}(3, 2)];
    c = [edges{idx}(1, 3); edges{idx}(2, 3); edges{idx}(3, 3)];
    d = [edges{idx}(1, 4); edges{idx}(2, 4); edges{idx}(3, 4)];

    pfunc = @(t) a(3) .* t.^3 + b(3) .* t.^2 + c(3) .* t + d(3);
    xfunc = @(t) a(1) .* pfunc(t).^3 + b(1) .* pfunc(t).^2 + c(1) .* pfunc(t) + d(1);
    yfunc = @(t) a(2) .* pfunc(t).^3 + b(2) .* pfunc(t).^2 + c(2) .* pfunc(t) + d(2);

    fplot(xfunc, yfunc, [0 T]);
end

hold off;

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