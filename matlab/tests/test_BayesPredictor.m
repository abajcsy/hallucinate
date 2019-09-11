%% This file tests the FULLY BAYESIAN predictor
clf
clear all

% Variance on Gaussian observation model.
sigma1 = pi/8;
sigma2 = pi/8;

% Known goal locations (in m). 
goals = {[1, 1], [1, -1]}; %{[1, tan(sigma1)], [1, -tan(sigma2)]};

% Grid structure definitions
gridMin = [-2,-2];          % Lower & upper bounds of grid (in m)
gridMax = [2,2];
gridDims = [41,41];         % Num grid cells in X and Y dimension. 

% Set the prior over goal 1 and goal 2.
prior = [0.5, 0.5];

% Create the predictor. 
predictor = BayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, gridDims);
%predictor.Pu_given_x_g(2, [21,21], 2);

% --- debugging ------ %
%predictor.plot_Pu_given_x();
% -------------------- %

% Initial human state (in m).
x0 = [0,0];
v = 0.6;

% Prediction horizon. 
gString = createGrid(gridMin, gridMax, gridDims);
dt = gString.dx(1)/v;
T = 2;                                  % horizon in (seconds)
H = floor(T/dt);                % horizon in (timesteps)

% Predict!
preds = predictor.predict(x0, H, v, dt, gString);

% Normalize preds
predsNorm = renormalizeProbs(preds, x0, v, dt, gString, H);

% Plot
figure(2)
set(gcf,'color','w');

epsilon = 0.0;
for t=1:H+1
    p = preds{t};
    pNorm = predsNorm{t};
    sumsump = sum(sum(p))
    sumsumpNorm = sum(sum(pNorm))
    
    % For visualizing all distributions
    figure(1)
    pcolor(gString.xs{2}, gString.xs{1}, p);
    colorbar
    caxis([0,1])
    title('original');
    
    figure(2)
    pcolor(gString.xs{2}, gString.xs{1}, pNorm);
    colorbar
    caxis([0,1])
    title('normalized');
    
    % (Unnormalized) For visualizing thresholded. 
    figure(3)
    p1 = (p > epsilon);
    [~,c1] = contour(gString.xs{2}, gString.xs{1}, p1, [0,0.1]);
    c1.LineWidth = 2;
    c1.Color = 'b';
    grid on
    
    hold on
    v = 0.6;
    r = v*(t-1)*dt + 0.05; 
    h = rectangle('Position',[x0(1)-r x0(2)-r 2*r 2*r],'Curvature',[1,1], 'EdgeColor', 'r');
    title(strcat('original t=', num2str((t-1)*dt)));
    
    % (Normalized) thresholded
    figure(4)
    p2 = (pNorm > epsilon);
    [M,c2] = contour(gString.xs{2}, gString.xs{1}, p2, [0,0.1]);
    c2.LineWidth = 2;
    c2.Color = 'b';
    grid on
    
    hold on
    v = 0.6;
    r = v*(t-1)*dt + 0.05; 
    h2 = rectangle('Position',[x0(1)-r x0(2)-r 2*r 2*r],'Curvature',[1,1], 'EdgeColor', 'r');
    title(strcat('normalized t=', num2str((t-1)*dt)));
    
    delete(c1);
    delete(c2);
    delete(h);
    delete(h2);
end

%% Sicc: renormalize them preds. Make sure you are getting rid of 
% the extra states you cant reach
function normPreds = renormalizeProbs(preds, x0, v, dt, g, H)
    normPreds = preds;
    for t = 2:H+1
        p = normPreds{t};
        
        % Compute FRS stuff.
        r = v*(t-1)*dt;
        frs = (g.xs{1} - x0(1)).^2 + (g.xs{2} - x0(2)).^2 - r^2;
        validFRSIndicies = find(frs <= 0);
        validPIndicies = find(p > 0);
        invalidIndicies = setdiff(validPIndicies,validFRSIndicies);
        
        % Sum the probabilities of all invalid indicies to redistribute
        % pmass.
        probRedistribute = sum(p(invalidIndicies));
        
        normPreds{t}(invalidIndicies) = 0.0;
        normPreds{t}(validFRSIndicies) = normPreds{t}(validFRSIndicies) + ...
                                probRedistribute/length(validFRSIndicies);
    end
end