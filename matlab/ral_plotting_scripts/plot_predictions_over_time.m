clear all
clf

%% Load data.
repo = what('hallucinate');
folder = '/ral_data/';
filename = 'bayesian_pg10.5.mat';

load(strcat(repo.path, folder, filename));


%% Setup important variables.
totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, -2], [2, 2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

%% Color setup.
start_color = [132, 224, 223]/255.;
end_color = [0, 158, 155]/255.;
red = linspace(start_color(1), end_color(1), totalTime);
green = linspace(start_color(2), end_color(2), totalTime);
blue = linspace(start_color(3), end_color(3), totalTime);
alpha = linspace(0.3, 0.9, totalTime);
hcolor = linspace(0.8, 0, totalTime);

%% Thresholds.
% Bayesian Predictions: Percent of states to omit from visualization.
% Reachability Predictions: Percent of controls to omit. 
uThresh = 0.03; 

incr = 8;
for t=1:incr:totalTime
    
    figure(1);
    hold on;

    % Plot the predictions.
    preds = all_preds{t};
    xcurr = human_states{t};
    
    predIncr = 2;
    predColorR = linspace(hcolor(t), red(t), length(preds));
    predColorG = linspace(hcolor(t), green(t), length(preds));
    predColorB = linspace(hcolor(t), blue(t), length(preds));
    for pt=1:predIncr:length(preds)
        curr_pred_color = [predColorR(pt),  predColorG(pt), predColorB(pt)];
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
            badIndicies = find(Yinner >2);
            Yinner(badIndicies) = [];
            Xinner(badIndicies) = [];
        end
        c = fill(Xinner,Yinner, curr_pred_color, ...
            'FaceAlpha', alpha(t), 'EdgeColor', 'none');
    end

    % Plot state of human.
    scatter(xcurr(1), xcurr(2), 60, [hcolor(t), hcolor(t), hcolor(t)], 'filled');
    
    % Plot goals.
    figure(1);
    scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
    scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');
    
    % Label goal 2.
    %text_pg1 = strcat("$P^t(g_1) = ", num2str(prior(1)), "$");
    %tpg1 = text(goals{1}(1)-0.4, goals{1}(2)-0.3, text_pg1, 'Color', 'r', 'Interpreter', 'Latex');
    tg1 = text(goals{1}(1)+0.1, goals{1}(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
    tg1.FontSize = 18;
    %tpg1.FontSize = 14;
    
    % Label goal 2.
    %text_pg2 = strcat("$P^t(g_2) = ", num2str(prior(2)), "$");
    %tpg2 = text(goals{2}(1)-0.4, goals{2}(2)+0.3, text_pg2, 'Color', 'r', 'Interpreter', 'Latex');
    tg2 = text(goals{2}(1)+0.1, goals{2}(2), '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
    tg2.FontSize = 18;
    %tpg2.FontSize = 14;

    xlim([-1, 4]);
    ylim([-2.5, 2.5]);
    set(gcf,'Position',[100 100 700 700]);
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

end

% Plot time annotations.
tt1 = text(human_states{1}(1)-0.4, human_states{1}(2), ...
        strcat('$t=',num2str(times(1),2),"$"), ...
        'Color', [hcolor(5), hcolor(5), hcolor(5)], 'Interpreter', 'Latex');
ttend = text(human_states{t}(1)+0.2, human_states{t}(2), ...
        strcat('$t=',num2str(times(t),2),"$"), ...
        'Color', [hcolor(t), hcolor(t), hcolor(t)], 'Interpreter', 'Latex');


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


