xDisc = 0.1;
yDisc = 0.1;
thetaDisc = pi / 45;
% vDisc = 0.05;
vDisc = 0.05;
tDisc = 0.5;

heurWeight = 1.25;

planner = LatticePlanner(xDisc, yDisc, thetaDisc, vDisc, tDisc, heurWeight);

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
planner.staticObsMap = OccupancyGrid(xDisc, yDisc, tDisc, ...
                                     xBounds(1), xBounds(2), ...
                                     yBounds(1), yBounds(2), ...
                                     tBounds(1), tBounds(2));
planner.dynObsMap = OccupancyGrid(xDisc, yDisc, tDisc, ...
                                  xBounds(1), xBounds(2), ...
                                  yBounds(1), yBounds(2), ...
                                  tBounds(1), tBounds(2));

startState = LatticePlannerState(0, 0, 0, 1, 0);
[succs, edges, costs] = planner.expand(startState);
fprintf("expanded %d successors\n", length(succs));

figure;
hold on;

for idx = 1:length(succs)
    a = [edges{idx}(1, 1); edges{idx}(2, 1); edges{idx}(3, 1)];
    b = [edges{idx}(1, 2); edges{idx}(2, 2); edges{idx}(3, 2)];
    c = [edges{idx}(1, 3); edges{idx}(2, 3); edges{idx}(3, 3)];
    d = [edges{idx}(1, 4); edges{idx}(2, 4); edges{idx}(3, 4)];

    pfunc = @(t) a(3) .* t.^3 + b(3) .* t.^2 + c(3) .* t + d(3);
    xfunc = @(t) a(1) .* pfunc(t).^3 + b(1) .* pfunc(t).^2 + c(1) .* pfunc(t) + d(1);
    yfunc = @(t) a(2) .* pfunc(t).^3 + b(2) .* pfunc(t).^2 + c(2) .* pfunc(t) + d(2);

    fplot(xfunc, yfunc, [0 tDisc]);
end

hold off;
