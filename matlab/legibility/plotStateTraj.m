%% Plots the states that result from optimal control sequence.
function plotStateTraj(uopt, human, times, tminIdx, goals, trueGoalIdx, grid_min, grid_max)
    figure(2);
    hold on

    % Keep track of all states so we can draw a line plot.
    X = zeros(1,length(times(1:tminIdx)));
    Y = zeros(1,length(times(1:tminIdx)));
    PG = zeros(1,length(times(1:tminIdx)));

    % Setup colors.
    startColor = [79., 0., 128.]/255.;
    endColor = [255., 143., 255.]/255.;
    r = linspace(startColor(1), endColor(1), length(times(1:tminIdx)));
    g = linspace(startColor(2), endColor(2), length(times(1:tminIdx)));
    b = linspace(startColor(3), endColor(3), length(times(1:tminIdx)));
    xcurr = human.x;

    % Record state.
    X(1) = xcurr(1);
    Y(1) = xcurr(2);
    PG(1) = xcurr(3);

    % Plot first point.
    color = [r(1), g(1), b(1)];
    plot3(xcurr(1), xcurr(2), xcurr(3), '-o', 'color', color, 'markeredgecolor', color, 'markerfacecolor', color);
    for t=2:length(times(1:tminIdx))
        % Apply control to dynamics.
        human.updateState(uopt{t}, times(t), human.x);
        xcurr = human.x;

        % record state.
        X(t) = xcurr(1);
        Y(t) = xcurr(2);
        PG(t) = xcurr(3);

        % Plot point and connection between pts.
        color = [r(t), g(t), b(t)];
        p = plot3([X(t-1), X(t)], [Y(t-1), Y(t)], [PG(t-1), PG(t)], '-o', ...
                'Color', color, ...
                'markeredgecolor', color, ...
                'markerfacecolor', color);
        p.LineWidth = 2;
        
        % add timestamps
        txt = num2str(times(t));
        t = text(X(t)+0.05, Y(t)+0.05, PG(t)+0.05, txt);
        t.Color = color;
    end

    % Plot goals (red is ground truth, grey is other goal).
    plot3(goals{trueGoalIdx}(1), goals{trueGoalIdx}(2), 0.5, '-o', ...
                'Color', 'r', ...
                'markeredgecolor', 'r', ...
                'markerfacecolor', 'r');
    g1Txt = strcat('g', num2str(trueGoalIdx));
    t1 = text(goals{trueGoalIdx}(1), goals{trueGoalIdx}(2), 0.55, g1Txt);
    t1.FontSize = 12;
    t1.Color = 'r';

    if trueGoalIdx == 2
        otherGoalIdx = 1; 
    else
        otherGoalIdx = 2;
    end
    plot3(goals{otherGoalIdx}(1), goals{otherGoalIdx}(2), 0.5, '-o', ...
                'Color', 'k', ...
                'markeredgecolor', 'k', ...
                'markerfacecolor', 'k');
            
    g2Txt = strcat('g', num2str(otherGoalIdx));
    t2 = text(goals{otherGoalIdx}(1), goals{otherGoalIdx}(2), 0.55, g2Txt);
    t2.FontSize = 12;
    t2.Color = 'k';
    
    grid on;
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    zlim([grid_min(3), grid_max(3)]);
    
    xlabel('x');
    ylabel('y');
    zlabel('P(g = g1)');
end