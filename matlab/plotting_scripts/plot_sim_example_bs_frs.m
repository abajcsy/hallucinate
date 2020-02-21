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
figure(1)
hold on
plotGoals(humanGoals, humanGoalColor);
plotGoalTol(robotGoals, robotGoalTol, robotGoalTolColor);
plotGoals(robotGoals, robotGoalColor);
scatter(2, 0, 250, 'r', 'LineWidth', 2);

%% Plot the predictions.
for i=1:length(human_pos)
    preds = squeeze(occupancy_maps(i,:,:,:));
    x = human_pos(i,:);
    
    if i == 1 %|| i == 10 %i == 7 || i == 15
        plotPredictions(x, preds, numPredsToSkip, predsColor, predsAlpha(i), g);
    end
end

%% Plot the human trajectory.
for i=1:1 %length(human_pos)
    x = human_pos(i,:);
    if i >= 1 && i < 3
        plotHuman(x, [0,0,0]);
    else
        plotHuman(x, [humanColor(i), humanColor(i), humanColor(i)]);
    end
    %plotHuman(x, [humanColor(i), humanColor(i), humanColor(i)]);
end

%% Plot the robot trajectory.
cIdx = 1;
for i=1:numStatesToSkip:length(position)
    if i > minStep && i < maxStep
        plotRobot(position(i,:), heading(i,:), ...
            [0, 64, 4]/255.);
    else
        plotRobot(position(i,:), heading(i,:), ...
            [0, robotColor(cIdx), 0]);
    end
	%plotRobot(position(i,:), heading(i,:), ...
    %    [robotColor(cIdx), robotColor(cIdx), robotColor(cIdx)]);
    cIdx = cIdx + 1;
end


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
        h{i} = scatter(g(1), g(2), 150, 'filled', ...
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
                'MarkerFaceColor', color, 'MaxHeadSize', 6.0, ...
                'ShowArrowHead', 'on', 'AutoScaleFactor', 0.5, ...
                'LineWidth', 2); 

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
	h = scatter(x(1), x(2), 150, 'filled', ...
        'MarkerEdgeColor', color,...
        'MarkerFaceColor', color); 
end

%% Plot human predictions.
function h = plotPredictions(x, preds, numPredsToSkip, color, alpha, g)
    idx = 1;
    for i=4:numPredsToSkip:length(preds)
        p = preds(:,:,i);
        
        [c, h] = contour(g.xs{1}, g.xs{2}, p, [0,1]);
        h.EdgeColor = color(idx,:);
        h.LineWidth = 1;
        X = c(1,2:end);
        Y = c(2,2:end);
        %delete(h);
        h = fill(X,Y,color(idx,:),'FaceAlpha', alpha , 'EdgeColor', 'none');

%         [~, h] = contourf(g.xs{1}, g.xs{2}, p, [0,1]);
%         colormap([1,1,1; color(idx,:)]);
        
%         for xs=1:length(p)
%             for ys=1:length(p)
%                 if p(xs, ys) > 0
%                     scatter(xs, ys, 5, color(idx,:), 'filled');
%                 end
%             end
%         end
        
        
%         [~, h] = contour(g.xs{1}, g.xs{2}, p, [0,0.01]);
%         h.EdgeColor = color(idx,:);
%         h.LineWidth = 1.5;
         idx = idx + 1;
    end
    
end
