% Get the parameters of the reachability problem.
%params = dubinsCarGaussianHuman();

% Create the predictor and make first prediction.
%predictor = HumanPredictor(params);

% Grab the predictions (value functions) and times.
%[preds, times] = predictor.getPredictions();xz

clc
clear all

% Create occupancy grid.
gmin = [0,0];
gmax = [2,2];
N = [3,2];
g = createGrid(gmin, gmax, N);
valueFun = shapeSphere(g, [1; 1.5], 1);
occuGrid = OccupancyGrid(N(1), N(2), 1, ...
    gmin(1),gmax(1),gmin(2),gmax(2));

% Convert to 2D binary occupancy map representation.
occuGrid.fromValueFuns(g, valueFun, 1, 0);
