params = scenario1_splines();

% FRS loading parameters.
% Assuming starting in the root directory of the project.
pathToFRSDir = './matlab/frs';
priorProb = 0;
horizon = 2;

predGrid2D = proj(params.predGrid, params.predGrid.xs, [0,0,1]);
dynObsMap = OccupancyGrid(params.xDisc, params.yDisc, params.dt, ...
    params.lowEnv(1), params.upEnv(1), ...
    params.lowEnv(2), params.upEnv(2), ...
    params.tMin, params.tMax, predGrid2D, pathToFRSDir);

% Load from file.
if ~dynObsMap.loadFromFile(priorProb, horizon)
    fprintf('Error: failed to load FRS!\n');
end

% Draw the resulting predictions.
figure;
hold on;

xlim([params.lowEnv(1), params.upEnv(1)]);
ylim([params.lowEnv(2), params.upEnv(2)]);

xH = [0; 0];
vH = [0.5; 0.25];

f = scatter([], []);
for t = dynObsMap.tMin:dynObsMap.tDisc:dynObsMap.tMax
    fprintf('At t = %f\n', t);
    
    xs = [];
    ys = [];
    
    % Simulate the human moving at a constant velocity.
    dynObsMap.setHumanState(xH);
    scatter([xH(1)], [xH(2)], 'xr');
    xH = xH + dynObsMap.tDisc * vH;
    
    for x = params.lowEnv(1):dynObsMap.xDisc:params.upEnv(1)
        for y = params.lowEnv(2):dynObsMap.yDisc:params.upEnv(2)
            % Get the occupancy at this (x, y).
            if dynObsMap.getData(x, y, t) > 0
                xs = [xs, x];
                ys = [ys, y];
            end
        end
    end   
    
    f.set('Xdata', xs);
    f.set('Ydata', ys);
    drawnow;
    pause(1);
end


hold off;
