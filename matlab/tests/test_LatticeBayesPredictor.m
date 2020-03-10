%% This file tests the FULLY BAYESIAN predictor
clf
clear all

% Variance on Gaussian observation model.
sigma1 = pi/4;
sigma2 = pi/4;

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; %{[1, tan(sigma1)], [1, -tan(sigma2)]};

% Grid structure definitions
gridMin = [-4,-4];          % Lower & upper bounds of grid (in m)
gridMax = [4,4];

% Set the prior over goal 1 and goal 2.
prior = [0.9, 0.1];

% Grid cell size.
r = 0.1;

% Create the predictor. 
predictor = LatticeBayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, r);

% Initial human state (in m).
x0 = [0,0];
v = 0.6;
dt = r / v;

% Prediction horizon. 
% gString = createGrid(gridMin, gridMax, gridDims);
% dt = gString.dx(1)/v;
T = 3;                                  % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)

% Predict!
preds = predictor.predict(x0, H);

% eps = 0.0000000001;
eps = 0.001;

figure(2);
hold on
for t=1:H+1
    xs = [];
    ys = [];
    ps = [];    
    p = preds{t};
    for s = predictor.states
        ss = s{1};
        [x, y] = predictor.simToReal(ss);
        xs = [xs, x];
        ys = [ys, y];
        ps = [ps, p(ss(1), ss(2))];
%         if p(ss(1), ss(2)) > eps
%             ps = [ps, 1];
%             fprintf('(%d, %d) --> (%f, %f)\n', ...
%                     ss(1), ss(2), x, y);
%         else
%             ps = [ps, 0];
%         end
    end
    
    sum(ps);
    
%     scatter(xs, ys, 30 * ones(1, length(ys)), 1 - ps, 'filled', 'MarkerEdgeColor', 'k');
%     colormap('gray');

%     sz = 30 * ones(1, length(ys));
%     %sz = 10;
% 	scatter(xs, ys, sz, ps, 'filled', 'MarkerEdgeColor', 'none');
%     scatter(goals{1}(1), goals{1}(2), 'r', 'filled');
%     scatter(goals{2}(1), goals{2}(2), 'r', 'filled');
%     xlim([gridMin(1), gridMax(1)]);
%     ylim([gridMin(2), gridMax(2)]);
%     colormap('gray');
%     colorbar
%     caxis([0 max(ps)]);

    titleString = strcat('Static param, t=', num2str(t*dt), 's');
    title(titleString);

    [X, Y] = predictor.getLatticeMeshgrid();
    
    P = zeros(size(X));
    for i = 1:predictor.rows
        for j = 1:predictor.cols
            P(i, j) = 1*(p(i, j) > eps) + 0*(p(i, j) <= eps);
        end
    end
    
    [M, c] = contour(X, Y, P, [1, 1]);
    c.LineWidth = 2;
    c.EdgeColor = [(195-t*10)/255, (191-t*10)/255, (191-t*10)/255]; 
    grid on
    pause(0.1);
end


% % Normalize preds
% predsNorm = renormalizeProbs(preds, x0, v, dt, gString, H);
% 
% % Plot
% figure(2)
% set(gcf,'color','w');
% 
% epsilon = 0.0;
% for t=1:H+1
%     p = preds{t};
%     pNorm = predsNorm{t};
%     sumsump = sum(sum(p))
%     sumsumpNorm = sum(sum(pNorm))
%     
%     % For visualizing all distributions
%     figure(1)
%     pcolor(gString.xs{2}, gString.xs{1}, p);
%     colorbar
%     caxis([0,1])
%     title('original');
%     
%     figure(2)
%     pcolor(gString.xs{2}, gString.xs{1}, pNorm);
%     colorbar
%     caxis([0,1])
%     title('normalized');
%     
%     % (Unnormalized) For visualizing thresholded. 
%     figure(3)
%     p1 = (p > epsilon);
%     [~,c1] = contour(gString.xs{2}, gString.xs{1}, p1, [0,0.1]);
%     c1.LineWidth = 2;
%     c1.Color = 'b';
%     grid on
%     
%     hold on
%     v = 0.6;
%     r = v*(t-1)*dt + 0.05; 
%     h = rectangle('Position',[x0(1)-r x0(2)-r 2*r 2*r],'Curvature',[1,1], 'EdgeColor', 'r');
%     title(strcat('original t=', num2str((t-1)*dt)));
%     
%     % (Normalized) thresholded
%     figure(4)
%     p2 = (pNorm > epsilon);
%     [M,c2] = contour(gString.xs{2}, gString.xs{1}, p2, [0,0.1]);
%     c2.LineWidth = 2;
%     c2.Color = 'b';
%     grid on
%     
%     hold on
%     v = 0.6;
%     r = v*(t-1)*dt + 0.05; 
%     h2 = rectangle('Position',[x0(1)-r x0(2)-r 2*r 2*r],'Curvature',[1,1], 'EdgeColor', 'r');
%     title(strcat('normalized t=', num2str((t-1)*dt)));
%     
%     delete(c1);
%     delete(c2);
%     delete(h);
%     delete(h2);
% end