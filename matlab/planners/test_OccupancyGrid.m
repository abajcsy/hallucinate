% Occupancy grid discretization.
xDisc = 0.15;
yDisc = 0.15;
tDisc = 0.5;

% Bounds on x coordinates.
xMin = 0;
xMax = 2;

% Bounds on y coordinates.
yMin = 0;
yMax = 2;

% Bounds on the time.
tMin = 0;
tMax = 1;

grid = OccupancyGrid(xDisc, yDisc, tDisc, xMin, xMax, yMin, yMax, tMin, ...
                     tMax);
grid.times = [0, 0.5, 1];

for x = 0.75:0.15:1.5
    for y = 0.5:0.15:1.4
        grid.setData(x, y, 0, 1);
    end
end

[kB, kA] = grid.closestIndexToTime(0)
[kB, kA] = grid.closestIndexToTime(0.1)
[kB, kA] = grid.closestIndexToTime(0.5)
[kB, kA] = grid.closestIndexToTime(0.55)
[kB, kA] = grid.closestIndexToTime(1.)
[kB, kA] = grid.closestIndexToTime(1.3)