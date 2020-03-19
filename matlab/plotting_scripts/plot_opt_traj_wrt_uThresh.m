clf
figure(4);
hold on

%% Setup directories for data.
%path = "/home/abajcsy/workspaces/hybrid_ws/src/hallucinate/data_uopt_x00/";
path = "/home/abajcsy/workspaces/hybrid_ws/src/hallucinate/data_uopt_x24/";

dirs = ["uopt_0thresh.mat", "uopt_001thresh.mat", "uopt_01thresh.mat", ...
        "uopt_05thresh.mat", "uopt_09thresh.mat"];
uthreshs = [0,0.01,0.1,0.5,0.9];

%% Setup colors.
numTrajs = length(dirs);

% startColor = [79., 0., 128.]/255.;
% endColor = [255., 143., 255.]/255.;
% r = linspace(startColor(1), endColor(1), numTrajs);
% g = linspace(startColor(2), endColor(2), numTrajs);
% b = linspace(startColor(3), endColor(3), numTrajs);

r = [0, 122, 188, 239, 255]/255.;
g = [63, 81, 80, 86, 166]/255.;
b = [92, 149, 144, 117, 0]/255.;

%% Plot all the trajectories.
for i=1:numTrajs
    color = [r(i), g(i), b(i)];
    load(strcat(path, dirs(i)));
    plotStateTraj(uopt, human, times, tminIdx, dt, color, uthreshs(i), i-1);
end

%% Plot goals (red is ground truth, grey is other goal).
rectangle('Position',[goals{trueGoalIdx}(1)-goalSetRad ...
                      goals{trueGoalIdx}(2)-goalSetRad ...
                      goalSetRad*2 ...
                      goalSetRad*2],...
                      'Curvature',1, ...
                      'FaceColor',[1, 0.67, 0.67],...
                      'EdgeColor',[1, 0.67, 0.67],...
                      'LineWidth',1);
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
rectangle('Position',[goals{otherGoalIdx}(1)-goalSetRad ...
                  goals{otherGoalIdx}(2)-goalSetRad ...
                  goalSetRad*2 ...
                  goalSetRad*2],...
                  'Curvature',1, ...
                  'FaceColor',[0.67, 0.67, 0.67],...
                  'EdgeColor',[0.67, 0.67, 0.67],...
                  'LineWidth',1);
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


%% Plots the states that result from optimal control sequence.
function plotStateTraj(uopt, human, times, tminIdx, dt, color, uthresh, idx)

    % Keep track of all states so we can draw a line plot.
    X = zeros(1,length(times(1:tminIdx)));
    Y = zeros(1,length(times(1:tminIdx)));
    PG = zeros(1,length(times(1:tminIdx)));

    xcurr = human.x;
    
    % marker size.
    sz = 3;

    % Record state.
    X(1) = xcurr(1);
    Y(1) = xcurr(2);
    PG(1) = xcurr(3);

    % Plot first point.
    plot3(xcurr(1), xcurr(2), xcurr(3), '-o', 'color', color, ...
        'markeredgecolor', color, 'markerfacecolor', color, 'markersize', sz);
    
    % Add annotation for threshold 
    txt = strcat('\delta=', num2str(uthresh));
    tp = text(-3.5 + 1.5*idx, -3, PG(1)+0.05, txt);
    tp.Color = color;
    tp.FontSize = 11;
    
    for t=2:length(times(1:tminIdx))
        % Apply control to dynamics.
        human.updateState(uopt{t}, dt, human.x);
        xcurr = human.x;

        % record state.
        X(t) = xcurr(1);
        Y(t) = xcurr(2);
        PG(t) = xcurr(3);

        % Plot point and connection between pts.
        p = plot3([X(t-1), X(t)], [Y(t-1), Y(t)], [PG(t-1), PG(t)], '-o', ...
                'Color', color, ...
                'markeredgecolor', color, ...
                'markerfacecolor', color, ...
                'markersize', sz);
        p.LineWidth = 1;
    end
    
    % Add annotation for final prob.
    txt2 = strcat('P(g1)=', num2str(PG(end), '%0.2f'));
    tp2 = text(-4 + 1.7*idx, -3.5, PG(end)+0.05, txt2);
    tp2.Color = color;
    tp2.FontSize = 8;
    
    % Add annotation for final time.
    txt2 = strcat('t_{min}=', num2str(times(end), '%0.1f'));
    tp2 = text(-3.8 + 1.7*idx, -4, PG(end)+0.05, txt2);
    tp2.Color = color;
    tp2.FontSize = 8;

%     %add timestamps
%     txt = strcat('t=', num2str(times(t)));%, ', p=', num2str(PG(t)));
%     %txt = strcat('e=', num2str(uthresh));
%     tp = text(X(t)-0.5, Y(t)-0.5, PG(t)+0.05, txt);
%     tp.Color = color;
%     tp.FontSize = 10;
end