%% This file tests the DYNAMIC-GOAL FULLY BAYESIAN predictor
clf
clear all

% Variance on Gaussian observation model.
sigma1 = pi/4;
sigma2 = pi/4;

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; %{[1, tan(sigma1)], [1, -tan(sigma2)]};

% Grid structure definitions
gridMin = [-4,-4];          % Lower & upper bounds of grid (in m)
gridMax = [4,4];

% Set the prior over goal 1 and goal 2.
prior = [0.9, 0.1];

% Grid cell size.
r = 0.1;

% HMM model parameter.
hmmParam = 0.6;

% Create the predictor. 
predictor = DynamicLatticeBayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, r, hmmParam);

% Initial human state (in m).
x0 = [0,0];
v = 0.6;
dt = r / v;

% Prediction horizon. 
% gString = createGrid(gridMin, gridMax, gridDims);
% dt = gString.dx(1)/v;
T = 3;                          % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)

% Predict!
preds = predictor.predict(x0, H);

% eps = 0.0000000001;
eps = 0.001;

figure(2);
hold on

for t=1:H+1
    xs = [];
    ys = [];
    ps = [];    
    p = preds{t};
    for s = predictor.states
        ss = s{1};
        [x, y] = predictor.simToReal(ss);
        xs = [xs, x];
        ys = [ys, y];
        ps = [ps, p(ss(1), ss(2))];
    end
    
    sum(ps)
    
%     sz = 30 * ones(1, length(ys));
%     %sz = 10;
% 	scatter(xs, ys, sz, ps, 'filled', 'MarkerEdgeColor', 'none');
%     scatter(goals{1}(1), goals{1}(2), 'r', 'filled');
%     scatter(goals{2}(1), goals{2}(2), 'r', 'filled');
%     xlim([gridMin(1), gridMax(1)]);
%     ylim([gridMin(2), gridMax(2)]);
%     colormap('pink');
%     colorbar
%     caxis([0 max(ps)]);

    titleString = strcat('HMM alpha=', num2str(hmmParam), ', t=', num2str(t*dt), 's');
    title(titleString);
    
    [X, Y] = predictor.getLatticeMeshgrid();
    
    P = zeros(size(X));
    for i = 1:predictor.rows
        for j = 1:predictor.cols
            P(i, j) = 1*(p(i, j) > eps) + 0*(p(i, j) <= eps);
        end
    end
    
    [M, c] = contour(X, Y, P, [1, 1]);
    c.LineWidth = 2;
    c.EdgeColor = [1, (191-t*10)/255, (195-t*10)/255]; 
    grid on
    pause(0.1);
end