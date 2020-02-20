clear all
clf
%close all

%% Load data.
path = "/home/abajcsy/hybrid_ws/src/hallucinate/data_for_paper/sim_example_bafrs/";
%filename = "pos_heading.mat";
filename = "all_data.mat";
load(strcat(path, filename)); % human_pos, occupancy_maps, vehicle_pos_heading

% Grab position and heading information about robot.
position = vehicle_pos_heading(:, 1:2);
heading = vehicle_pos_heading(:, 3);
human_pos = human_pos(1:10, :);

%% Plotting stuff.

% World bounds.
xlimits = [-2.5,4];
ylimits = [-2.5,2.5];
numStatesToSkip = 15;
numPredsToSkip = 20;

numRobotColors = ceil(length(position)/numStatesToSkip);
numHumanColors = length(human_pos);
numPredsColors = ceil(length(occupancy_maps)/numPredsToSkip);

% Known human goal locations. 
humanGoals = {[2,2], [2,-2]};  
robotGoals = {[-2,0]};  
robotGoalTol = 0.3;

% Indicies of the states that represent the planning horizon
minStep = 1;
maxStep = 80;

%% Colors.
robotColor = linspace(0.4, 0.7, numRobotColors);
humanColor = linspace(0.5, 0.8, numHumanColors);
predsColor = [linspace(0, 0.9412, numPredsColors); ...
              linspace(0, 0.1333, numPredsColors); ...
              linspace(0, 0.6706, numPredsColors)].';
predsAlpha = linspace(0.15, 0.6, length(human_pos));
humanGoalColor = 'r';
robotGoalColor = [0, 191, 6]/255.;
robotGoalTolColor = [181, 245, 183]/255.;

% Create grid for plotting in real-world.
lowEnv = [-4;-4];
upEnv = [4;4];
N = [81; 81];
g = createGrid(lowEnv, upEnv, N);

%% Plot human and robot goals.
hold on
figure(3)
plotGoals(humanGoals, humanGoalColor);
plotGoalTol(robotGoals, robotGoalTol, robotGoalTolColor);
plotGoals(robotGoals, robotGoalColor);
scatter(2, 0, 150, 'r');

%% Plot the human trajectory.
for i=1:length(human_pos)
    x = human_pos(i,:);
    if i >= 1 && i < 3
        plotHuman(x, [0,0,0]);
    else
        plotHuman(x, [humanColor(i), humanColor(i), humanColor(i)]);
    end
    %plotHuman(x, [humanColor(i), humanColor(i), humanColor(i)]);
end

%% Plot the robot
plotRobot(position(1,:), heading(1,:), 'k');

xlim(xlimits);
ylim(ylimits);
xticks([]);
yticks([]);
box on
set(gcf,'color','w');

%% Plot goals
function h = plotGoals(goals, color)
    h = {};
    for i=1:length(goals)
        g = goals{i};
        h{i} = scatter(g(1), g(2), 60, 'filled', ...
            'MarkerEdgeColor', color, ...
            'MarkerFaceColor', color);
    end
end

%% Plot goal tolernace. 
function h = plotGoalTol(goals, goalTol, color)
    g = goals{1};
    h = rectangle('Position', [g(1)-goalTol, g(2)-goalTol, goalTol*2, goalTol*2], ...
            'Curvature', 1, ...
            'FaceColor', color, ...
            'EdgeColor', color);
end

%% Plot robot.
function h = plotRobot(x, theta, color)
    h = quiver(x(1), x(2), cos(theta), sin(theta), 'o', 'Color', color, ...
                'MarkerSize', 7, 'MarkerEdgeColor', color, ...
                'MarkerFaceColor', color, 'MaxHeadSize', 5.0, ...
                'ShowArrowHead', 'on', 'AutoScaleFactor', 0.4, ...
                'LineWidth', 1); 

%     h = {};
% 	h{1} = scatter(x(1), x(2), 'filled', ...
%         'MarkerEdgeColor', color,...
%         'MarkerFaceColor', color); 
%     
%     % Rotation matrix.
%     R = [cos(theta) -sin(theta); 
%          sin(theta) cos(theta)];
%      
%     % Heading pt.
%     hpt = [0.5; 0];
%     hptRot = R*hpt + x;
%     h{2} = plot([x(1) hptRot(1)], [x(2) hptRot(2)], 'Color', color, 'LineWidth', 1.5);
end

%% Plot human. 
function h = plotHuman(x, color)
	h = scatter(x(1), x(2), 'filled', ...
        'MarkerEdgeColor', color,...
        'MarkerFaceColor', color); 
end

