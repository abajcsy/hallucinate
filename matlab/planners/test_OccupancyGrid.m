% Occupancy grid discretization.
xDisc = 0.15;
yDisc = 0.15;

% Bounds on x coordinates.
xMin = 0;
xMax = 2;

% Bounds on y coordinates.
yMin = 0;
yMax = 2;

grid = OccupancyGrid(xDisc, yDisc, xMin, xMax, yMin, yMax);

for x = 0.75:0.15:1.5
    for y = 0.5:0.15:1.4
        grid.setData(x, y, 1);
    end
end
