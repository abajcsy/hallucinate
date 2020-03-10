clear all;
close all;
clc;

% Figure parameters
start_color = [0,0,0]/255.;
end_color = [240, 34, 171]/255.; % Magenta
xLimits = [-2 2.5];
yLimits = [-2.5 2.5];
zlimits = [0, 1];
widthMeters = 8;
heightMeters = 8;
MarkerSize = 20;
ticksize = 14;
fontsize = 30;
alpha_min = 0.3;
alpha_max = 0.8;

% Parameters to plot     
priors = [0.8];
delta = [0.1]; %[0.3];
% indices_to_plot = [20, 40, 60];
indices_to_plot = [25]; %[20];
subindices_to_plot = {[2, 10, 19]};
num_plots = length(indices_to_plot);

% Create a temporary grid for the 2D sets
grid_min = [-4; -4; zlimits(1)];  % Lower corner of computation domain
grid_max = [4; 4; zlimits(1) + 0.05];     % Upper corner of computation domain
N = [81; 81; 41];           % Number of grid points per dimension
gtemp = createGrid(grid_min, grid_max, N);

%% Start generating the plots
 % Load the reachability based predictions
path = '/home/abajcsy/hybrid_ws/src/hallucinate';
filename = strcat(path, '/data_for_paper/example2/', 'predictions_reahcability_prior_', num2str(priors(1)), '_delta_', num2str(delta(1)), '.mat');
load(filename) % Fetch 'g', 'data'

% Start Plotting 
for i=1:num_plots
  % Total number of sets to plot
  total_plots = 1 + length(subindices_to_plot{i}) + 1; % Start set, end set and the subindices sets
  
  % Create the figure
  f = figure('Position', 2*[100, 100, widthMeters*100*0.6, heightMeters*100*0.6]); 
  set(gcf,'Color','w') 
  hold on

  % Plot the initial condition in 3D 
  %plot3(0, 0, priors(1), 'o', 'MarkerSize', MarkerSize, 'MarkerfaceColor', start_color);
  initData = shapeSphere(g, [0;0;priors(1)], 0.1);
  visSetIm(g, initData, start_color, 0.0);
  
  % Plot the initial condition in 2D 
  gtemp2D = createGrid(grid_min(1:2), grid_max(1:2), N(1:2));
  initData2D = shapeSphere(gtemp2D, [0;0], 0.1);
  [~, c] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, initData2D, [0,0.01]);
  c.LineWidth = 2.5;
  c.EdgeColor = start_color;
  %plot3(0, 0, 0, 'MarkerSize', MarkerSize, 'MarkerfaceColor', start_color);
  
  % Plot the goal 1 in 2D.
  g1Data2D = shapeSphere(gtemp2D, [2;2], 0.1);
  [~, c2] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, g1Data2D, [0,0.01]);
  c2.LineWidth = 2.5;
  c2.EdgeColor = 'r';
  
  % Plot the goal 1 in 2D.
  g2Data2D = shapeSphere(gtemp2D, [2;-2], 0.1);
  [~, c3] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, g2Data2D, [0,0.01]);
  c3.LineWidth = 2.5;
  c3.EdgeColor = 'r';
  
  % Plot the FRS for subindices
  alphas = linspace(alpha_min, alpha_max, length(subindices_to_plot{i})+1);
  for j=1:length(subindices_to_plot{i})
    % Color of the set
    current_color = (start_color * (1 - ((j+1)/total_plots))) +  (end_color * ((j+1)/total_plots));
    
    % Plot the set in 3D
    h{j} = visSetIm(g, data(:, :, :, subindices_to_plot{i}(j)), current_color, 0.0);
    h{j}.FaceAlpha = alphas(j);
    
    % Plot the set in 2D
    [gtemp2D, data2D] = proj(gtemp, data(:, :, :, subindices_to_plot{i}(j)), [0,0,1], 'min');
    [~, l{j}] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, data2D, [0,0.01]);
    l{j}.LineWidth = 2.5;
    l{j}.EdgeColor = current_color;
    %l{j} = visSetIm(gtemp, data(:, :, :, subindices_to_plot{i}(j)), current_color, 0.0);
    %l{j}.FaceAlpha = alphas(j);
    
  end
  
  % Plot the end set in 3D
  h{j+1} = visSetIm(g, data(:, :, :, indices_to_plot(i)), end_color, 0.0);
  h{j+1}.FaceAlpha = alphas(j+1);
  
  % Plot the end set in 2D
  [gtemp2D, data2D] = proj(gtemp, data(:, :, :, indices_to_plot(i)), [0,0,1], 'min');
  [~, l{j+1}] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, data2D, [0,0.01]);
  l{j+1}.LineWidth = 2.5;
  l{j+1}.EdgeColor = end_color;
  %l{j+1} = visSetIm(gtemp, data(:, :, :, indices_to_plot(i)), end_color, 0.0);
  %l{j+1}.FaceAlpha = alphas(j+1);
    
  c = camlight;
  c.Position = [119.5625 119.2780 15.3565];
  
  % Set the view angle
  az = -7.6;
  el = 34.2;
  view(az,el)

  % Set the right figure size and limits
  xlim(xLimits)
  ylim(yLimits)
  zlim(zlimits)
  set(gca, 'xtick', [-4, -2, 0, 2, 4]);
  set(gca, 'ytick', [-4, -2, 0, 2, 4]);
  set(gca, 'ztick', [0, 0.25, 0.5, 0.75, 1.0]);
  set(gcf, 'Position',  0.5*[100, 100, 1500, 1200])
  set(gca, 'fontsize', ticksize);
  xlabel('$p_x$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
  ylabel('$p_y$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
  zlabel('$P(g_1)$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
  
  % Save the figure
  filename = strcat(path, '/data_for_paper/example2/', 'sets_prior_', num2str(priors(1)), '_delta_', num2str(delta(1)), '_index_', num2str(indices_to_plot(i)));
  saveas(f, strcat(filename, '.fig'))
  saveas(f, strcat(filename, '.png'))
  saveas(f, strcat(filename, '.pdf'))

end














% % Plot the reachability predictor
% s1 = contour(g2D.xs{1}, g2D.xs{2}, pred2D, [0, 0], 'color', colors{1}, 'LineWidth', LineWidth,'LineStyle',linestyle{1});
% 
% % Plot the Bayesian predictor
% s2 = contour(X, Y, P, [1, 1], 'color', colors{2}, 'LineWidth', LineWidth,'LineStyle',linestyle{2});
% 
% % Plot the full FRS
% R = v * time_s;
% pos = [- R, - R, 2*R, 2*R];
% s3 = rectangle('Position',pos,'Curvature',1.0,'LineStyle',linestyle{3},'LineWidth',LineWidth, 'EdgeColor', colors{3});
% 
% % Plot the goals
% R = goal_radius;
% pos = [goals{1}(1)- R, goals{1}(2)- R, 2*R, 2*R];
% s4 = rectangle('Position',pos,'Curvature',1.0,'FaceColor',colors{4},'LineStyle','none');
% pos = [goals{2}(1)- R, goals{2}(2)- R, 2*R, 2*R];
% s5 = rectangle('Position',pos,'Curvature',1.0,'FaceColor',colors{4},'LineStyle','none');
% 
% % Set the right figure size and limits
% xlim(xLimits);
% ylim(yLimits);
% set(gcf, 'Position',  0.5*[100, 100, widthMeters*100*0.6, heightMeters*100*0.6])
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% box on


