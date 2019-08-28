%% Clear old figure plotting and variables.
clf 
clc 
clear all

%% Load the experimental setup.
params = dubinsCarGaussianHuman();

%% Create human predictor.
predictor = HumanPredictor(params);

%% Variables that we are saving out.
timestamps = {};
humanStates = {};
humanPreds = {};
predsTimes = {};

%% Simulation loop.
hold on

% Create the first prediction.
predictor.updatePredictions();
[preds, times] = predictor.getPredictions();

% project to 2D
data2D = [];
for i=1:length(times)
    [~, projPred] = proj(params.predGrid, preds(:,:,:,i), [0,0,1], 'min');
    data2D = cat(3, data2D, projPred);
end

% Get the initial state of the simulated human.cvf
xHnext = params.xH0(1:2);
xHgoal = [1.5; -0.5];
hh = [];

for t=0:params.T-1
    % If human is close enough to goal position, terminate. 
    if norm(xHnext - xHgoal) <= 0.2
        break;
    end
    
    % Record data.
    timestamps{end+1} = t;
    humanStates{end+1} = xHnext;
    humanPreds{end+1} = data2D;
    predsTimes{end+1} = times;

    % ----- plotting ------ %
%     if ~isempty(hh)
%         delete(hh{1});
%     end
%     hh = plotAgent(xHnext, 'b');        
%     xlim([params.lowEnv(1),params.upEnv(1)]);
%     ylim([params.lowEnv(2),params.upEnv(2)]);
%     pause(0.1);
    % --------------------- %

    % Get most recent measurement of where the person is and what action
    % they applied.
    [xHnext, uHcurr] = params.simHuman.simulateAction(params.simDt);

    % Prediction step.
    predictor.updateState(xHnext, uHcurr);
    predictor.updatePredictions();
    [preds, times] = predictor.getPredictions();
    
    % Project and save. 
    data2D = [];
    for i=1:length(times)
        [~, projPred] = proj(params.predGrid, preds(:,:,:,i), [0,0,1], 'min');
        data2D = cat(3, data2D, projPred);
    end
end

% Save out data.
save('introspective_p09.mat', 'timestamps', 'humanStates', 'humanPreds', 'predsTimes');

%% Plots human or robot.
% Inputs:
%   x [vector]  - 3D/4D state of agent
% Ouput:
%   c   - handle for figure
function c = plotAgent(x, color)
    c = {};
    c{1} = plot(x(1), x(2), 'ko','MarkerSize', 8, ...
        'MarkerEdgeColor', color, 'MarkerFaceColor', color);

    % Plot heading.
    center = x(1:2);

    if length(x) >= 3
        % Rotation matrix.
        R = [cos(x(3)) -sin(x(3)); 
             sin(x(3)) cos(x(3))];
        % Heading pt.
        hpt = [0.2; 0];
        hptRot = R*hpt + center;
        p2 = plot([center(1) hptRot(1)], [center(2) hptRot(2)], color, 'LineWidth', 1.5);
        p2.Color(4) = 1.0;
        c{2} = p2;
    end
end