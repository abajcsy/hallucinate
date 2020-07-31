clear all


%% Load data.
repo = what('hallucinate');
folder = '/ral_data/';

filename = 'bayesian_modelledg_static_pg10.5_delta0.01.mat';
%filename = 'bayesian_unknowng_static_pg10.5_delta0.01.mat';
%filename = 'bayesian_subopt_static_pg10.5_delta0.01.mat';
pred_method = 'bayes';

%filename = 'reach_modelledg_static_pg10.5_gamma0.5_delta0.01.mat';
%filename = 'reach_unknowng_static_pg10.5_gamma0.5_delta0.01.mat';
%filename = 'reach_subopt_static_pg10.5_gamma0.5_delta0.01.mat';

%filename = 'reach_modelledg_static_pg10.5_gamma0.4_delta0.02.mat';
%filename = 'reach_unknowng_static_pg10.5_gamma0.4_delta0.02.mat';
%filename = 'reach_subopt_static_pg10.5_gamma0.4_delta0.02.mat';
%pred_method = 'reach';

%pred_method = 'frs';

load(strcat(repo.path, folder, filename));

%% Setup important variables.
totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

incr = 6; % time increments in which to plot everything.

%% Plot the posterior over time?
figNum = 1;
plot_posterior(times, posteriors, figNum, incr, totalTime);

%% Color setup.
bayes_start_color = [132, 224, 223]/255.;
bayes_end_color = [0, 158, 155]/255.;
reach_start_color = [224, 132, 198]/255.;
reach_end_color = [158, 0, 129]/255.;
frs_start_color = [224, 224, 224]/255.;
frs_end_color = [138, 138, 138]/255.;

if strcmp(pred_method, 'bayes')
    start_color = bayes_start_color;
    end_color = bayes_end_color;
elseif strcmp(pred_method, 'reach')
    start_color = reach_start_color;
    end_color = reach_end_color;
elseif strcmp(pred_method, 'frs')
    start_color = frs_start_color;
    end_color = frs_end_color;
else
    error('Illegal prediction method!');
end

red = linspace(start_color(1), end_color(1), totalTime);
green = linspace(start_color(2), end_color(2), totalTime);
blue = linspace(start_color(3), end_color(3), totalTime);
alpha = linspace(0.3, 0.5, totalTime);
hcolor = linspace(0.7, 0, totalTime);

%% Thresholds.
% Bayesian Predictions: Percent of states to omit from visualization.
% Reachability Predictions: Percent of controls to omit. 
%uThresh = 0.03; 
uThresh = 0.02;
%uThresh = 0.01; 

figure(1);
hold on;

% plot unmodelled obstacle.
obs = rectangle('Position', [-0.5,-1.5,1.5,1.5], 'Curvature', 1);
obs.LineStyle = '-.';
obs.EdgeColor = [0.6,0.6,0.6];
obs.LineWidth = 1;
obs.FaceColor = [0.9,0.9,0.9];

for t=1:incr:totalTime

    % Plot the predictions.
    preds = all_preds{t};
    xcurr = human_states{t};
    
    if strcmp(pred_method, 'bayes')
        finalTime = length(preds);
    else
        finalTime = length(preds(1,1,1,:));
    end
    
    predIncr = 2;
    predColorR = linspace(hcolor(t), red(t), finalTime);
    predColorG = linspace(hcolor(t), green(t), finalTime);
    predColorB = linspace(hcolor(t), blue(t), finalTime);
    
    for pt=1:predIncr:finalTime
        curr_pred_color = [predColorR(pt),  predColorG(pt), predColorB(pt)];
        
        if strcmp(pred_method, 'bayes')
            curr_pred = preds{pt};
            [opt_eps, P, X, Y] = ...
                compute_likely_states(curr_pred, predictor, uThresh);

            % Plot prediction contour.
            [M, c] = contour(X, Y, P, [1, 1]);
            c.LineWidth = 1;
            c.EdgeColor = curr_pred_color;
            Xinner = M(1,2:end);
            Yinner = M(2,2:end);
            if pt > 0
                badIndicies = find(Yinner > 2);
                Yinner(badIndicies) = [];
                Xinner(badIndicies) = [];
            end
            
            if pt == 9 && t == 25
                badIndicies = [26:1:33];
                Xinner(badIndicies) = [];
                Yinner(badIndicies) = [];
            end
            
            c = fill(Xinner,Yinner, curr_pred_color, ...
                'FaceAlpha', alpha(t), 'EdgeColor', 'none');
            %patch(Xinner,Yinner, curr_pred_color);
        elseif strcmp(pred_method, 'reach')
            % Plot prediction contour.
            dimsToRemove = [0 0 1];
            [g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
            [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
            c.LineWidth = 1;
            c.EdgeColor = curr_pred_color;
            Xinner = M(1,2:end);
            Yinner = M(2,2:end);
            c = fill(Xinner,Yinner, curr_pred_color, ...
                'FaceAlpha', alpha(t), 'EdgeColor', 'none');
        else
            v = 0.6;
            frs_rad = v*all_pred_times{t}(pt);
            pos = [xcurr(1)-frs_rad, xcurr(2)-frs_rad, frs_rad*2, frs_rad*2];
            full_frs = rectangle('Position', pos, 'Curvature', 1);
            full_frs.LineStyle = '-';
            full_frs.EdgeColor = curr_pred_color;
            full_frs.LineWidth = 1;
            full_frs.FaceColor = [curr_pred_color, alpha(t)];
        end
    end

    % Plot state of human.
    scatter(xcurr(1), xcurr(2), 60, [hcolor(t), hcolor(t), hcolor(t)], 'filled');
    
    % Plot goals.
    figure(1);
    scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
    scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');
    
    % Label goal 1.
    %text_pg1 = strcat("$P^t(g_1) = ", num2str(prior(1)), "$");
    %tpg1 = text(goals{1}(1)-0.4, goals{1}(2)-0.3, text_pg1, 'Color', 'r', 'Interpreter', 'Latex');
    tg1 = text(goals{1}(1)+0.3, goals{1}(2)-0.1, '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
    tg1.FontSize = 25;
    %tpg1.FontSize = 14;
    
    % Label goal 2.
    %text_pg2 = strcat("$P^t(g_2) = ", num2str(prior(2)), "$");
    %tpg2 = text(goals{2}(1)-0.4, goals{2}(2)+0.3, text_pg2, 'Color', 'r', 'Interpreter', 'Latex');
    tg2 = text(goals{2}(1)+0.3, goals{2}(2)+0.1, '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
    tg2.FontSize = 25;
    %tpg2.FontSize = 14;

%     xlim([-1, 4]);
%     ylim([-2.5, 2.5]);
%     set(gcf,'Position',[100 100 700 700]);

    xlim([-2, 4]);
    ylim([-2.5, 2.5]);
    set(gcf,'Position',[100 100 700 600]);
    set(gcf,'color','w');
    whitebg('w');
    
    box on
    set(gca,'xcolor','k','ycolor','k', ...
        'xtick',[], 'xticklabel',[], ...
        'ytick',[], 'yticklabel',[])
    set(gcf,'color','w');
    
%     % Plot posterior. 
%     figure(3);
%     hold on
%     plot(times, posteriors, 'r-o', 'LineWidth', 3);
%     xlabel("$time$", 'Interpreter', 'Latex');
%     ylabel("$P(g_1)$", 'Interpreter', 'Latex');
%     ylim([0,1]);
%     grid on

%     tt = text(human_states{t}(1)+0.4, human_states{t}(2), ...
%         strcat('$t=',num2str(times(t),2),"$"), ...
%         'Color', [hcolor(t), hcolor(t), hcolor(t)], 'Interpreter', 'Latex');
end

% Plot time annotations.
tt1 = text(human_states{1}(1)-0.8, human_states{1}(2), ...
        strcat('$t=',num2str(times(1),2),"$"), ...
        'Color', [hcolor(5), hcolor(5), hcolor(5)], 'Interpreter', 'Latex');
tt1.FontSize = 18;

% ttend = text(human_states{t}(1)+0.8, human_states{t}(2), ...
%         strcat('$t=',num2str(times(t),2),"$"), ...
%         'Color', [hcolor(t), hcolor(t), hcolor(t)], 'Interpreter', 'Latex');
% ttend.FontSize = 12;

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

%% Plots the posterior over time.
function plot_posterior(times, posteriors, figNum, incr, totalTime)
    figure(figNum);
    
    hold on;
    times_sampled = zeros(1,floor(totalTime/incr));
    posteriors_sampled = zeros(1,floor(totalTime/incr));
    str_times = cell(1,floor(totalTime/incr));
    ii = 1;
    for t=1:incr:totalTime
        times_sampled(ii) = times(t);
        str_times{ii} = num2str(times_sampled(ii), 2);
        posteriors_sampled(ii) = posteriors(t);
        ii = ii+1;
    end
    pp = plot(times_sampled, posteriors_sampled, '-o');
    pp.LineWidth = 5;
    pp.MarkerFaceColor = 'r';
    pp.Color = 'r';
    
    xlim([0, times(end)]);
    ylim([0,1]);
    set(gcf,'Position',[100 100 400 100]);
    set(gcf,'color','w');
    whitebg('w');
    xlab = xlabel('$t (sec)$', 'Interpreter', 'Latex');
    ylab = ylabel('$b^t(g_1)$', 'Interpreter', 'Latex');
    xlab.FontSize = 20;
    ylab.FontSize = 20;
    
    box on
    set(gca,'xcolor','k','ycolor','k', ...
        'xtick',[], 'xticklabel', str_times, ...
        'ytick',[], 'yticklabel',[])
    set(gcf,'color','w');
end


