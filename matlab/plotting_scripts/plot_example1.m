clear all;
close all;
clc;

% Figure parameters
LineWidth = 1.5;
magenta = [240, 34, 171]/255.;
teal = [0, 153, 153]/255.;
babyBlue = [84, 167, 222]/255.;
grey = [158, 158, 158]/255.;
red = [255, 0, 0]/255;
darkRed = [153, 0, 0]/255.;
colors = {magenta, teal, grey, red, red}; % In the order for reachability, bayesian, full FRS predictors, and the goals.
linestyle = {'-', '-', '-'}; % In the order for reachability, bayesian, and full FRS predictors.
xLimits = [-3 3];
yLimits = [-3 3];
widthMeters = 8;
heightMeters = 8;
MarkerSize = 5;
goal_radius = 0.2;

% Parameters to plot     
priors = [0.5, 0.9];
delta = [0.05, 0.1, 0.2];

num_priors = length(priors);
num_eps = length(delta);

confidence = 0.05;

%% Start generating the plots
for i=1:num_priors
  % Load the Bayesian predictor
  path = '/home/abajcsy/hybrid_ws/src/hallucinate';
  filename = strcat(path, '/data_for_paper/example1/', 'predictions_bayesian_prior_', num2str(priors(i)), '.mat');
  load(filename) % Fetch 'preds', 'predictor', 'goals', 'v', 'time_s', 'dt'
  
  for j=1:num_eps
    % Load the reachability based predictions
    filename = strcat(path, '/data_for_paper/example1/', 'predictions_reahcability_prior_', num2str(priors(i)), '_delta_', num2str(delta(j)), '.mat');
    load(filename) % Fetch 'g2D', 'pred2D'
    
    % Compute an appropriate epsilon for the Bayesian predictor
    [eps, P, X, Y] = compute_the_right_eps(pred2D, predictor, preds, 'discard_unlikely_states', confidence); %delta(j));
    
    % Start Plotting 
    f = figure;
    hold on;
    
    % Plot the position of the human
    scatter(0, 0, 'k', 'filled');
    
    % Plot the full FRS
    R = v * time_s + 0.2;
    pos = [- R, - R, 2*R, 2*R];
    s3 = rectangle('Position',pos,'Curvature',1.0,'LineStyle',linestyle{3},'LineWidth',LineWidth, 'EdgeColor', colors{3});
    
    % Plot the reachability predictor
    s1 = contour(g2D.xs{1}, g2D.xs{2}, pred2D, [0, 0], 'color', colors{1}, 'LineWidth', LineWidth,'LineStyle',linestyle{1});
    
    % Plot the Bayesian predictor
    s2 = contour(X, Y, P, [1, 1], 'color', colors{2}, 'LineWidth', LineWidth,'LineStyle',linestyle{2});
    
    % Plot the goals
    R = goal_radius;
    pos = [goals{1}(1)- R, goals{1}(2)- R, 2*R, 2*R];
    s4 = rectangle('Position',pos,'Curvature',1.0,'FaceColor',colors{4},'LineStyle','none');
    pos = [goals{2}(1)- R, goals{2}(2)- R, 2*R, 2*R];
    s5 = rectangle('Position',pos,'Curvature',1.0,'FaceColor',colors{5},'LineStyle','none');
    
    % Set the right figure size and limits
    xlim(xLimits);
    ylim(yLimits);
    set(gcf, 'Position',  0.5*[100, 100, widthMeters*100*0.6, heightMeters*100*0.6])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    box on
    
    % Save the figure
    filename = strcat(path, '/data_for_paper/example1/', 'sets_prior_', num2str(priors(i)), '_delta_', num2str(delta(j)));
    saveas(f, strcat(filename, '.fig'))
    saveas(f, strcat(filename, '.png'))
    saveas(f, strcat(filename, '.pdf'))
    
  end
end
