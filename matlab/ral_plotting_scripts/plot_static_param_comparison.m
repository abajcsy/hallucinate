clear all
clf

%% ======= PLOT BAYESIAN PREDICTION ======= %

%% ------- Prior = 0.5 --------- %
% Load data.
repo = what('hallucinate');
folder = '/ral_data/';
filename = 'bayesian_static_pg10.5_1pred.mat';
load(strcat(repo.path, folder, filename));

totalTime = length(all_preds);

% Setup thresholding param.
uThresh = 0.03;

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

preds = all_preds{1};
xcurr = human_states{1};

% [bayes] Color setup.
bayes_start_color = [0, 77, 73]/255.;
bayes_middle_color = [15, 166, 146]/255.; 
bayes_end_color = [0, 199, 189]/255.; %[47, 219, 119]/255.;
b_red = [linspace(bayes_start_color(1), bayes_middle_color(1), length(preds)/2), ...
            linspace(bayes_middle_color(1), bayes_end_color(1), length(preds)/2+1)];
b_green = [linspace(bayes_start_color(2), bayes_middle_color(2), length(preds)/2), ...
            linspace(bayes_middle_color(2), bayes_end_color(2), length(preds)/2+1)];
b_blue = [linspace(bayes_start_color(3), bayes_middle_color(3), length(preds)/2), ...
            linspace(bayes_middle_color(3), bayes_end_color(3), length(preds)/2+1)];
alpha = linspace(0.6, 0.8,length(preds));

figure(1);
sfh1 = subplot(2,2,1);
hold on;

% Label title and yaxis.
sft1 = title('Intent-Driven Bayes', 'Interpreter', 'Latex');
sft1.FontSize = 14;
yt1 = ylabel('$b^0(\lambda = g_1) = 0.5$', 'Interpreter', 'Latex');
yt1.FontSize = 14;

% Plot the predictions.
predIncr = 3;
small_dt = 0.1667;
vel = 0.6;
for pt=1:predIncr:length(preds)
    curr_pred_color = [b_red(pt),  b_green(pt), b_blue(pt)];
    curr_pred = preds{pt};
    [opt_eps, P, X, Y] = ...
        compute_likely_states(curr_pred, predictor, uThresh);
    
    % Plot full forward reachable set (FRS) = dt*vel
    DT = (pt-1)*small_dt;
    frs_rad = DT*vel;
    frs_color = [alpha(pt), alpha(pt), alpha(pt)]; %curr_pred_color;
    frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
        'Curvature',1, 'EdgeColor', frs_color);
    frs.LineStyle = '--';

    % Plot prediction contour.
    [M, c] = contour(X, Y, P, [1, 1]);
    c.LineWidth = 2;
    c.EdgeColor = curr_pred_color;
end

% Plot state of human.
scatter(xcurr(1), xcurr(2), 60, 'k', 'filled');

% Plot goals.
figure(1);
scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');

% Label goal 2.
tg1 = text(goals{1}(1)+0.1, goals{1}(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.1, goals{2}(2), '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-2, 2.5]);
ylim([-2.5, 2.5]);
set(gcf,'Position',[100 100 700 700]);
set(gcf,'color','w');
whitebg('w');

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
set(gcf,'color','w');
hold off;

%% ------- Prior = 0.9 --------- %
clear all

% Setup thresholding param.
uThresh = 0.03;

% Load data.
repo = what('hallucinate');
folder = '/ral_data/';
filename = 'bayesian_static_pg10.9_1pred.mat';
load(strcat(repo.path, folder, filename));

totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

preds = all_preds{1};
xcurr = human_states{1};

% [bayes] Color setup.
bayes_start_color = [0, 77, 73]/255.;
bayes_middle_color = [15, 166, 146]/255.; 
bayes_end_color = [0, 199, 189]/255.; %[47, 219, 119]/255.;
b_red = [linspace(bayes_start_color(1), bayes_middle_color(1), length(preds)/2), ...
            linspace(bayes_middle_color(1), bayes_end_color(1), length(preds)/2+1)];
b_green = [linspace(bayes_start_color(2), bayes_middle_color(2), length(preds)/2), ...
            linspace(bayes_middle_color(2), bayes_end_color(2), length(preds)/2+1)];
b_blue = [linspace(bayes_start_color(3), bayes_middle_color(3), length(preds)/2), ...
            linspace(bayes_middle_color(3), bayes_end_color(3), length(preds)/2+1)];
alpha = linspace(0.6, 0.8,length(preds));

% Plot
sfh2 = subplot(2,2,3);
hold on;

% Label y-axis.
yt2 = ylabel('$b^0(\lambda = g_1) = 0.9$', 'Interpreter', 'Latex');
yt2.FontSize = 14;

% Plot the predictions.
predIncr = 3;
small_dt = 0.1667;
vel = 0.6;
for pt=1:predIncr:length(preds)
    curr_pred_color = [b_red(pt),  b_green(pt), b_blue(pt)];
    curr_pred = preds{pt};
    [opt_eps, P, X, Y] = ...
        compute_likely_states(curr_pred, predictor, uThresh);
    
    % Plot full forward reachable set (FRS) = dt*vel
    DT = (pt-1)*small_dt;
    frs_rad = DT*vel;
    frs_color = [alpha(pt), alpha(pt), alpha(pt)]; %curr_pred_color;
    frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
        'Curvature',1, 'EdgeColor', frs_color);
    frs.LineStyle = '--';

    % Plot prediction contour.
    [M, c] = contour(X, Y, P, [1, 1]);
    c.LineWidth = 2;
    c.EdgeColor = curr_pred_color;
end

% Plot state of human.
scatter(xcurr(1), xcurr(2), 60, 'k', 'filled');

% Plot goals.
figure(1);
scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');

% Label goal 2.
tg1 = text(goals{1}(1)+0.1, goals{1}(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.1, goals{2}(2), '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-2, 2.5]);
ylim([-2.5, 2.5]);
set(gcf,'Position',[100 100 700 700]);
set(gcf,'color','w');
whitebg('w');

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
set(gcf,'color','w');

% Setup the position of the second figure.
sfh2.Position = sfh2.Position + [0,0.1,0,0];

%% ======= PLOT REACH PREDICTION ======= %

%% ------- Prior = 0.5 --------- %

clear all

% Load data.
repo = what('hallucinate');
folder = '/ral_data/';
filename = 'reach_static_pg10.5_gamma0.5_delta0.03_1pred.mat';
load(strcat(repo.path, folder, filename));

totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

preds = all_preds{1};
xcurr = human_states{1};

% [reachability] Color setup.
reach_start_color = [97, 17, 90]/255.;
reach_end_color = [235, 19, 216]/255.;
r_red = linspace(reach_start_color(1), reach_end_color(1), length(preds(1,1,1,:)));
r_green = linspace(reach_start_color(2), reach_end_color(2), length(preds(1,1,1,:)));
r_blue = linspace(reach_start_color(3), reach_end_color(3), length(preds(1,1,1,:)));
alpha = linspace(0.8,0.6,length(preds(1,1,1,:)));

sfh3 = subplot(2,2,2);
hold on;

% Label title and yaxis.
sft3 = title('BA-FRS', 'Interpreter', 'Latex');
sft3.FontSize = 14;

dimsToRemove = [0 0 1];

% Plot the predictions.
predIncr = 3;
small_dt = 0.1667;
vel = 0.6;
for pt=1:predIncr:length(preds(1,1,1,:))
    curr_pred_color = [r_red(pt),  r_green(pt), r_blue(pt)];
   
    % Plot full forward reachable set (FRS) = dt*vel
    DT = (pt-1)*small_dt;
    frs_rad = DT*vel;
    frs_color = [alpha(pt), alpha(pt), alpha(pt)]; %curr_pred_color;
    frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
        'Curvature',1, 'EdgeColor', frs_color);
    frs.LineStyle = '--';

    % Plot prediction contour.
    [g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
    [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
    c.LineWidth = 2;
    c.EdgeColor = curr_pred_color;
end

% Plot state of human.
scatter(xcurr(1), xcurr(2), 60, 'k', 'filled');

% Plot goals.
figure(1);
scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');

% Label goal 2.
tg1 = text(goals{1}(1)+0.1, goals{1}(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.1, goals{2}(2), '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-2, 2.5]);
ylim([-2.5, 2.5]);
set(gcf,'Position',[100 100 700 700]);
set(gcf,'color','w');
whitebg('w');

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
set(gcf,'color','w');

%0.3347 0.3412
%sfh3.Position = sfh3.Position + [-0.1,-0.025,0,0];

%% ------- Prior = 0.9 --------- %
clear all

% Load data.
repo = what('hallucinate');
folder = '/ral_data/';
filename = 'reach_static_pg10.9_gamma0.5_delta0.03_1pred.mat';
load(strcat(repo.path, folder, filename));

totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

preds = all_preds{1};
xcurr = human_states{1};

% [reachability] Color setup.
reach_start_color = [97, 17, 90]/255.;
reach_end_color = [235, 19, 216]/255.;
r_red = linspace(reach_start_color(1), reach_end_color(1), length(preds(1,1,1,:)));
r_green = linspace(reach_start_color(2), reach_end_color(2), length(preds(1,1,1,:)));
r_blue = linspace(reach_start_color(3), reach_end_color(3), length(preds(1,1,1,:)));
alpha = linspace(0.8,0.6,length(preds(1,1,1,:)));

sfh4 = subplot(2,2,4);
hold on;

dimsToRemove = [0 0 1];

% Plot the predictions.
predIncr = 3;
small_dt = 0.1667;
vel = 0.6;
for pt=1:predIncr:length(preds(1,1,1,:))
    curr_pred_color = [r_red(pt),  r_green(pt), r_blue(pt)];
   
    % Plot full forward reachable set (FRS) = dt*vel
    DT = (pt-1)*small_dt;
    frs_rad = DT*vel;
    frs_color = [alpha(pt), alpha(pt), alpha(pt)]; %curr_pred_color;
    frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
        'Curvature',1, 'EdgeColor', frs_color);
    frs.LineStyle = '--';

    % Plot prediction contour.
    [g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
    [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
    c.LineWidth = 2;
    c.EdgeColor = curr_pred_color;
end

% Plot state of human.
scatter(xcurr(1), xcurr(2), 60, 'k', 'filled');

% Plot goals.
figure(1);
scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');

% Label goal 2.
tg1 = text(goals{1}(1)+0.1, goals{1}(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.1, goals{2}(2), '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-2, 2.5]);
ylim([-2.5, 2.5]);
set(gcf,'Position',[100 100 700 700]);
set(gcf,'color','w');
whitebg('w');

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
set(gcf,'color','w');

% Setup the position of the second figure.
sfh4.Position = sfh4.Position + [-0.1,0.1,0,0];

sfh3 = subplot(2,2,2);
sfh3.Position = sfh3.Position + [-0.1,0,0,0];
hold off

% %% Title the figure.
% xpos = sfh3.Position(1);
% ypos = sfh3.Position(2)+sfh3.Position(3);
% t1 = text(xpos, ypos, '\textbf{Static} $\lambda^t$');%, 'Interpreter', 'Latex');
% t1.FontSize = 15;

%% Save!
repo = what('hallucinate');
savefolder = '/ral_imgs/';
savefilename = strcat('static_param_comparison_delta', num2str(uThresh),'.png');
saveas(gcf, strcat(repo.path, savefolder, savefilename));

%% Grab all the likely-enough predicted states.
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

