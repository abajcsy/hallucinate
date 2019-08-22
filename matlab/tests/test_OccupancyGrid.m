% Get the parameters of the reachability problem.
params = dubinsCarGaussianHuman();

% Create the predictor and make first prediction.
predictor = HumanPredictor(params);

% Grab the predictions (value functions) and times.
[preds, times] = predictor.getPredictions();

% Create occupancy grid. 
tDisc = (params.tMax - params.tMin)/params.dt;
occuGrid = OccupancyGrid(params.predGrid.N(1), params.predGrid.N(2), tDisc, ...
    params.lowEnv(1), params.upEnv(1), params.lowEnv(2), params.upEnv(2));

% Convert to 2D binary occupancy map representation.
occuGrid.fromValueFuns(params.predGrid, preds, times, 0);