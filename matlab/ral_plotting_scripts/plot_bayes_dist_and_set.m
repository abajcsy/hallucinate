clear all
close all
clf

%% Load data.
repo = what('hallucinate');
folder = '/ral_data/';
savefilename = 'bayesian_static_pg10.5_1pred.mat';

load(strcat(repo.path, folder, savefilename));

%% Setup important variables.
totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, -2], [2, 2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

%% Time slice to grab
initialSlice = 1;
predTimeSlice = 9;

preds = all_preds{initialSlice};
xcurr = human_states{initialSlice};

% Visualization lims.
% xlims = [-1, 2.5];
% ylims = [-2.5, 2.5];
xlims = [-1, 1.5];
ylims = [-1.5, 1.5];

%% Plot full distribution.
f1h = figure(1);
sfh1 = subplot(1,3,1);
hold on
xs = [];
ys = [];
ps = [];    
plotting_pred = preds{predTimeSlice};
for s = predictor.states
    ss = s{1};
    [x, y] = predictor.simToReal(ss);
    sprob = plotting_pred(ss(1), ss(2));
    
    % Remove all values on the boundary of the plotting (for clearner vis).
    if x <= xlims(1)+0.1 || y == ylims(1) || x >= xlims(2)-0.1 || y == ylims(2) ...
            || y < ylims(1)+0.2 || y > ylims(2)-0.2 
        continue;
    end
    if sprob > 0
        xs = [xs, x];
        ys = [ys, y];
        ps = [ps, sprob];
    end
end

%% Subfigure info.
subfigW = 0.3;
subfigH = 0.8;
offsetX = 0.1;
offsetY = 0.1;
titleSz = 17;

%% Setup colors.
num_p_vals = 1000;

near_zero_color = [235, 255, 235]/255.; %[255, 255, 255]/255.;
start_color = [229, 255, 228]/255.;
middle1_color = [190, 255, 191]/255.;
middle2_color = [47, 219, 119]/255.;
end_color = [25, 176, 145]/255.;

red = [linspace(near_zero_color(1), start_color(1), 10), ...
        linspace(start_color(1), middle1_color(1), num_p_vals/3.), ...
        linspace(middle1_color(1), middle2_color(1), num_p_vals/3.), ...
        linspace(middle2_color(1), end_color(1), num_p_vals/3.)];
green = [linspace(near_zero_color(2), start_color(2), 10), ...
        linspace(start_color(2), middle1_color(2), num_p_vals/3.), ...
        linspace(middle1_color(2), middle2_color(2), num_p_vals/3.), ...
        linspace(middle2_color(2), end_color(2), num_p_vals/3.)];
blue = [linspace(near_zero_color(3), start_color(3), 10), ...
        linspace(start_color(3), middle1_color(3), num_p_vals/3.), ...
        linspace(middle1_color(3), middle2_color(3), num_p_vals/3.), ...
        linspace(middle2_color(3), end_color(3), num_p_vals/3.)];
bayes_colormap = horzcat(red', green', blue');

% Setup colormap and colorbar.
colormap(bayes_colormap);
cbh = colorbar;
cbh.Ticks = [min(ps), max(ps)/2., max(ps)];
cbh.TickLabels = {'0', num2str(max(ps)/2.,2), num2str(max(ps),2)};

%% Plot the prediction distribution.
% Plot distribution! 
szone = 50;
sz = szone * ones(1, length(ys));
scatter(xs, ys, sz, ps, 'filled', 'MarkerEdgeColor', 'none');

%% Plot goals.
% scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
% scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');
scatter(1,-1, 100, 'r', 'filled');
scatter(1,1, 100, 'r', 'filled');
g2_txt = text(1,-1-0.2, '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
g2_txt.FontSize = 14;
g1_txt = text(1,1+0.2, '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
g1_txt.FontSize = 14;

%% Plot state of human.
scatter(xcurr(1), xcurr(2), szone, 'k', 'filled');

%% Setup axes and title.
t1 = title('$P(x^{1.8}_H \mid x^0_H)$', 'Interpreter', 'Latex');
t1.FontSize = titleSz;
xlim(xlims);
ylim(ylims);

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
whitebg('w');
set(gcf,'color','w');
%hold off;

%% ---------------- THRESHOLDED SET 1 ---------------- %%
sfh2 = subplot(1,3,2);
hold on 
delta = 0.0;
[opt_eps, P, X, Y] = compute_likely_states(plotting_pred, predictor, delta);

% Plot full forward reachable set (FRS) = dt*vel
DT = (predTimeSlice-1)*0.1667;
frs_rad = DT*0.6;
frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
    'Curvature',1, 'EdgeColor', [0.7,0.7,0.7]);
frs.LineStyle = '--';

% Plot prediction contour.
[M, c] = contour(X, Y, P, [1, 1]);
c.LineWidth = 2;
c.EdgeColor = end_color ;

% Plot goals.
% scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
% scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');
scatter(1,-1, 100, 'r', 'filled');
scatter(1,1, 100, 'r', 'filled');
g2_txt = text(1,-1-0.2, '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
g2_txt.FontSize = 14;
g1_txt = text(1,1+0.2, '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
g1_txt.FontSize = 14;

% Plot state of human.
scatter(xcurr(1), xcurr(2), szone, 'k', 'filled');

% Setup axes and title.
title_text = strcat('$P(x^{1.8}_H \mid x^0_H) >', num2str(delta), '$');
t1 = title(title_text, 'Interpreter', 'Latex');
t1.FontSize = titleSz;
xlim(xlims);
ylim(ylims);

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
whitebg('w');
set(gcf,'color','w');
%hold off;

%% ---------------- THRESHOLDED SET 2 ---------------- %%
sfh3 = subplot(1,3,3);
hold on 
delta = 0.1;
[opt_eps, P, X, Y] = compute_likely_states(plotting_pred, predictor, delta);

% Plot full forward reachable set (FRS) = dt*vel
DT = (predTimeSlice-1)*0.1667;
frs_rad = DT*0.6;
frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
    'Curvature',1, 'EdgeColor', [0.7,0.7,0.7]);
frs.LineStyle = '--';

% Plot prediction contour.
[M, c] = contour(X, Y, P, [1, 1]);
c.LineWidth = 2;
c.EdgeColor = end_color ;

% Plot goals.
% scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
% scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');
scatter(1,-1, 100, 'r', 'filled');
scatter(1,1, 100, 'r', 'filled');
g2_txt = text(1,-1-0.2, '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
g2_txt.FontSize = 14;
g1_txt = text(1,1+0.2, '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
g1_txt.FontSize = 14;

% Plot state of human.
scatter(xcurr(1), xcurr(2), szone, 'k', 'filled');

% Setup axes and title.
title_text = strcat('$P(x^{1.8}_H \mid x^0_H) >', num2str(delta), '$');
t1 = title(title_text, 'Interpreter', 'Latex');
t1.FontSize = titleSz;
xlim(xlims);
ylim(ylims);

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
whitebg('w');
set(gcf,'color','w');

%% Setup all positions
% set(f1h,'OuterPosition',[0, 0, 1000, 400]);
% set(sfh1, 'Position', [offsetX offsetY subfigW subfigH]);
% set(sfh2, 'Position', [offsetX+subfigW+0.1 offsetY subfigW subfigH]);
% set(sfh3, 'Position', [offsetX+subfigW*2+0.2 offsetY subfigW subfigH]);

set(f1h,'OuterPosition',[0, 0, 900, 400]);
set(sfh1, 'Position', [0.1, 0.1, 0.2, 0.7]);
set(sfh3, 'Position', [0.1+0.4+0.125, 0.1, 0.2, 0.7]);
set(sfh2, 'Position', [0.1+0.2+0.09, 0.1, 0.2, 0.7]);

%% Save!
savefolder = '/ral_imgs/';
savefilename = 'bayes_dist_and_set.png';
saveas(gcf, strcat(repo.path, savefolder, savefilename));

%% =========== HELPER FUNCTIONS ========== %%
% Grab all the likely-enough predicted states.
function [opt_eps, P, X, Y] = compute_likely_states(preds, predictor, ...
    delta_reachability)
    
    % Grid for Bayesian prediction
    [X, Y] = predictor.getLatticeMeshgrid();
    
    valid_indices = find(preds > 0);
    valid_data = preds(valid_indices);
    sorted_valid_data = sort(valid_data, 'descend');
    eps_index = find(cumsum(sorted_valid_data) > (1 - delta_reachability), 1, 'first');
    
    if isempty(eps_index)
        % if we can't find likely enough states, then we should take max. 
        eps_index = find(max(cumsum(sorted_valid_data)), 1, 'first');
    end
    
    opt_eps = sorted_valid_data(eps_index);
    
    % Compute the optimal predictions
    P = zeros(size(X));
    for r = 1:predictor.rows
        for c = 1:predictor.cols
            linIdx = sub2ind([predictor.rows,predictor.cols],r,c);
            if delta_reachability > 0.0
                P(r, c) = 1*(preds(linIdx) >= opt_eps) + 0*(preds(linIdx) < opt_eps);
            else
                P(r, c) = 1*(preds(linIdx) > 0.0);
            end
        end
    end
end
