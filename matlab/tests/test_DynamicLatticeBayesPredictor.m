%% This file tests the DYNAMIC-GOAL FULLY BAYESIAN predictor
clf
clear all

% Variance on Gaussian observation model.
sigma1 = pi/4;
sigma2 = pi/4;

% Known goal locations (in m). 
goals = {[2, -2], [2, 2]}; 

% Lower & upper bounds of grid (in m).
gridMin = [-4,-4];          
gridMax = [4,4];

% Grid cell size.
r = 0.1;

% Set the prior over goal 1 and goal 2.
prior = [0.5, 0.5];

% Discrete distribution over two goals: \delta{1,2}.
Delta = [0.5, 0.5];

% Probability of staying in same goal (in HMM model).
hmmParam = 0.5;

% Create the predictor. 
predictor = DynamicLatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
    gridMin, gridMax, r, hmmParam, Delta);

% Initial human state (in m).
x0 = [0,0];
v = 0.6;
dt = r / v;

% Prediction horizon. 
% gString = createGrid(gridMin, gridMax, gridDims);
% dt = gString.dx(1)/v;
T = 1.8;                          % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)

% Predict!
tStart = tic;
preds = predictor.predict(x0, H);
tEnd = toc(tStart); 
fprintf(strcat("Dynamic param pred time: ", num2str(tEnd), " s\n"));

% Reachability threshold.
uThresh = 0.01;

figure(1);
hold on

% Color setup.
start_color = [0, 96, 112]/255.;
end_color = [0, 186, 217]/255.;
red = linspace(start_color(1), end_color(1), H+1);
green = linspace(start_color(2), end_color(2), H+1);
blue = linspace(start_color(3), end_color(3), H+1);

for t=1:H+1

    % Get the current predictions, and the sufficiently likley volume.
    curr_pred = preds{t};
    [opt_eps, P, X, Y] = compute_likely_states(curr_pred, predictor, uThresh);
    %fprintf(strcat("opt_eps: ", num2str(opt_eps), "\n"));
    
    % Setup title.
    figure(1);
    titleString = strcat('Dynamic Param. HMM alpha=', num2str(hmmParam), ', t=', num2str(t*dt), 's');
    title(titleString);
    
    figure(1);
    hold on;
    
    %Plot prediction contour.
    [M, c] = contour(X, Y, P, [1, 1]);
    c.LineWidth = 2;
    c.EdgeColor = [red(t), green(t), blue(t)];
  
%     xs = [];
%     ys = [];
%     ps = [];    
%     p = preds{t};
%     for s = predictor.states
%         ss = s{1};
%         [x, y] = predictor.simToReal(ss);
%         xs = [xs, x];
%         ys = [ys, y];
%         ps = [ps, p(ss(1), ss(2))];
%     end
%     sz = 100; %30 * ones(1, length(ys));
%     hold on;
% 	  scatter(xs, ys, sz, ps, 'filled', 'MarkerEdgeColor', 'none');
%     colormap('pink');
%     colorbar
%     caxis([0 max(ps)]);
    
    % Plot goals.
    figure(1);
    scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
    scatter(goals{2}(1), goals{2}(2), 100, 'b', 'filled');
    
    %xlim([-0.5,0.5]);
    %ylim([-0.5,0.5]);
        
    xlim([gridMin(1), gridMax(1)]);
    ylim([gridMin(2), gridMax(2)]);
    set(gcf,'Position',[100 100 700 700]);
    set(gcf,'color','w');
    whitebg('k');
    grid on
    pause(0.1);
end

%% Grab all the likely-enough predicted states.
function [opt_eps, P, X, Y] = compute_likely_states(preds, predictor, ...
    delta_reachability)
    
    % Grid for Bayesian prediction
    [X, Y] = predictor.getLatticeMeshgrid();
    
    valid_indices = find(preds > 0);
    valid_data = preds(valid_indices);
    sorted_valid_data = sort(valid_data, 'descend');
    eps_index = find(cumsum(sorted_valid_data) > (1 - delta_reachability), 1, 'first');
    opt_eps = sorted_valid_data(eps_index);

    % Compute the optimal predictions
    P = zeros(size(X));
    for r = 1:predictor.rows
        for c = 1:predictor.cols
            linIdx = sub2ind([predictor.rows,predictor.cols],r,c);
            P(r, c) = 1*(preds(linIdx) >= opt_eps) + 0*(preds(linIdx) < opt_eps);
        end
    end
end