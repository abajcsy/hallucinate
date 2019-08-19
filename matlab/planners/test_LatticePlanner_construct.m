% xDisc = 0.25;
% yDisc = 0.25;
xDisc = 0.15;
yDisc = 0.15;
thetaDisc = pi / 45;
vDisc = 0.5;
tDisc = 0.25;

planner = LatticePlanner(xDisc, yDisc, thetaDisc, vDisc, tDisc);

% Set up the state bounds.
planner.stateBounds('x') = [-0.1, 3];
planner.stateBounds('y') = [-0.1, 3];
planner.stateBounds('v') = [0.1, 3];
planner.stateBounds('theta') = [-pi, pi];

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
goalXY = [1.5; 2];
goalTol = 0.2;
startStateCont = [0; 0; 0; 0.5; 0];
traj = planner.plan(startStateCont, goalXY, goalTol);

figure;
hold on;
axis equal;

drawTriangle([startStateCont(1); startStateCont(2)], startStateCont(3), 0.1);
scatter([goalXY(1)], [goalXY(2)]);
drawCircle(goalXY, goalTol);

for idx = 1:length(traj)-1
    firstState = traj{idx};
    secondState = traj{idx+1};

    drawTriangle([planner.discToCont(secondState.x, 1); ...
                  planner.discToCont(secondState.y, 2)], ...
                 planner.discToCont(secondState.theta, 3), ...
                 0.1);

    [a, b, c, d] = unicycleThirdOrderSpline(planner.discToCont(firstState.x, 1), ...
                                            planner.discToCont(firstState.y, 2), ...
                                            planner.discToCont(firstState.theta, 3), ...
                                            planner.discToCont(firstState.v, 4), ...
                                            planner.discToCont(secondState.x, 1), ...
                                            planner.discToCont(secondState.y, 2), ...
                                            planner.discToCont(secondState.theta, 3), ...
                                            planner.discToCont(secondState.v, 4), ...
                                            planner.discToCont(1, 5));

    xfunc = @(t) a(1) .* t.^3 + b(1) .* t.^2 + c(1) .* t + d(1);
    yfunc = @(t) a(2) .* t.^3 + b(2) .* t.^2 + c(2) .* t + d(2);

    vxfunc = @(t) 3 .* a(1) .* t.^2 + 2 .* b(1) .* t + c(1);
    vyfunc = @(t) 3 .* a(2) .* t.^2 + 2 .* b(2) .* t + c(2);

    fplot(xfunc, yfunc, [0 tDisc]);
end
hold off;

function tri = drawTriangle(origin, rot, sideLength)
    aLocal = [2*sideLength; 0; 1];
    bLocal = [0; -sideLength / sqrt(2); 1];
    cLocal = [0; sideLength / sqrt(2); 1];

    T = [cos(rot), -sin(rot), origin(1);
         sin(rot), cos(rot), origin(2);
         0, 0, 1];
    a = T * aLocal;
    b = T * bLocal;
    c = T * cLocal;

    tri = fill([a(1), b(1), c(1)], ...
               [a(2), b(2), c(2)], ...
               'b');
    set(tri, 'facealpha', 0.5);
end

function circ = drawCircle(origin, radius)
    th = 0:pi/50:2*pi;
    xunit = radius * cos(th) + origin(1);
    yunit = radius * sin(th) + origin(2);
    circ = plot(xunit, yunit);
end
