clf
clear all

% nine actions.
actions = [-pi, -(3*pi)/4, -pi/2, -pi/4, 0, pi/4, pi/2, (3*pi)/4, pi];
actionsStr = {'-pi', '-(3*pi)/4', '-pi/2', '-pi/4', '0', 'pi/4', 'pi/2', '(3*pi)/4', 'pi'};
%action_bounds = {[(7*pi)/8, -(7*pi)/8], ...
%                  [-(7*pi)/8; -(5*pi)/8], ...
%                  [-(5*pi)/8; -(3*pi)/8], ...
%                  [-(3*pi)/8; -(pi)/8], ...
%                  [-(pi)/8; (pi)/8], ...
%                  [(pi)/8; (3*pi)/8], ...
%                  [(3*pi)/8; (5*pi)/8], ...
%                  [(5*pi)/8; (7*pi)/8], ...
%                  [(7*pi)/8, -(7*pi)/8]};

action_bounds = {[-(7*pi)/8, (7*pi)/8], ...
                 [-(7*pi)/8; -(5*pi)/8], ...
                 [-(5*pi)/8; -(3*pi)/8], ...
                 [-(3*pi)/8; -(pi)/8], ...
                 [-(pi)/8; (pi)/8], ...
                 [(pi)/8; (3*pi)/8], ...
                 [(3*pi)/8; (5*pi)/8], ...
                 [(5*pi)/8; (7*pi)/8], ...
                 [-(7*pi)/8, (7*pi)/8]};


% standard dev of gaussian.
sigma = pi/8;

% goal locations (meters).
g1 = [1; -1];
g2 = [1; 1];

% state of human.
x = [0;0];

% 50/50 prior on each goal.
priorg1 = 0.5;

% action probabilities
action_probs_g1 = zeros(1,length(actions));
action_probs_g2 = zeros(1,length(actions));

% posteriors
posteriors = zeros(1,length(actions));

% comput probability for goal 1.
for i=1:length(actions)
    u = actions(i);
    ubounds = action_bounds{i};
    
    % compute P(u|x,g1) and P(u|x,g2)
    pug1 = PuGivenGoal(x, u, g1, sigma, ubounds, action_bounds);
    pug2 = PuGivenGoal(x, u, g2, sigma, ubounds, action_bounds);
    
    % compute P(g1|u,x) \propto P(u|x,g1)*P(g1)
    posteriorpug1 = (pug1 * priorg1)/(pug1 * priorg1 + pug2 *(1-priorg1));
    
    % printing..
    fprintf(strcat("P(",actionsStr{i},")=", num2str(pug1), "\n"));
    fprintf(strcat("   P(g1 | u, x) = ", num2str(posteriorpug1), "\n"));
    
    % plotting info...
    action_probs_g1(i) = pug1;
    action_probs_g2(i) = pug2;
    posteriors(i) = posteriorpug1;
end

% Plot action probabilities: P(u|x,g1)
figure(1)
plot_circle(actions, action_probs_g1, g1, g2);
title("P(u|x,g1)");

% Plot action probabilities: P(u|x,g2)
figure(2)
plot_circle(actions, action_probs_g2, g1, g2);
title("P(u|x,g2)");

% Plot posterior: P(g1|u,x)
figure(3)
hold on
plot(actions, posteriors, 'mo-')
plot([-pi, pi], [priorg1, priorg1], 'k--');
xticks(actions)
xticklabels(actionsStr)
set(gcf,'color','w');
grid on
title("Posterior for G1")
hold off

function pu = PuGivenGoal(x, u, goal, sigma, ubounds, action_bounds)
    
    % comput optimal action.
    uopt = atan2(goal(2)- x(2), goal(1) - x(1)); 
    
      %% ------- OPTION 1 ------ %
%     dToGoal = wrapToPi(u - uopt);
%     %dToGoal = u - uopt;
%    
%     c0 = 1/(sqrt(2*pi)*sigma);
%     uG1Diff = -(dToGoal)^2;
%     pu = c0 * exp(uG1Diff / (2*sigma^2));
%     % ------- OPTION 1 ------ %
    
%     %% ------- OPTION 2 ------ %
%     pd = makedist('Normal','mu',0,'sigma',sigma);
%     truncpd = truncate(pd, -pi, pi);
% 
%     % Find the integration bounds. 
%     bound1 = wrapToPi(uopt - ubounds(1));
%     bound2 = wrapToPi(uopt - ubounds(2));
% 
%     % Integrate on bounds. 
%     p = cdf(truncpd, [bound1, bound2]);
%     pu = abs(p(2) - p(1));
%     % ------- OPTION 2 ------ %
    
    %% ------- OPTION 2.5 (CORRECT VERSION) ------ %
    diff = abs(angdiff(u, uopt));
    
    % Find the integration bounds. 
    [trueBounds, mult] = getTrueBounds(diff, action_bounds);
    
    pd = makedist('Normal','mu',0,'sigma',sigma);
    truncpd = truncate(pd, -pi, pi);
    
    % Integrate on bounds. 
    p = cdf(truncpd, [trueBounds(1), trueBounds(2)]);
    pu = abs(p(2) - p(1))*mult;
    % ------- OPTION 2.5 (CORRECT VERSION) ------ %

%     %% ------- OPTION 3 ------ %
%     pd = makedist('Normal','mu',0,'sigma',sigma);
%     
%     % Find the integration bounds. 
%     bound1 = wrapToPi(uopt - ubounds(1));
%     bound2 = wrapToPi(uopt - ubounds(2));
%     
%     % Compute normalizer.
%     cdf_u_low = cdf(pd, -pi);
%     cdf_u_up = cdf(pd, pi);
%     trunc_normalizer = cdf_u_up - cdf_u_low;
%     
%     % Integrate on bounds. 
%     plow = cdf(pd, bound1);
%     pup = cdf(pd, bound2);
%     pu = abs(pup - plow)/trunc_normalizer;
%     
%     % Wrap prrobability.
%     if (abs(pup - plow) > pi)
%         pu = 1 - pu;
%     end
%     % ------- OPTION 3 ------ %
end

function [bounds, mult] = getTrueBounds(diff, action_bounds)
    mult = 1;
    if (0 <= diff && diff < pi / 8) 
        bounds = [0, pi/8];
        mult = 2;
    elseif (pi / 8 <= diff && diff < (3 * pi) / 8) 
        bounds = [pi/8, (3*pi)/8];
    elseif ((3 * pi) / 8 <= diff && diff < (5 * pi) / 8) 
        bounds = [(3*pi)/8, (5*pi)/8];
    elseif ((5 * pi) / 8 <= diff && diff < (7 * pi) / 8) 
        bounds = [(5*pi)/8, (7*pi)/8];
    elseif ((7 * pi) / 8 <= diff) 
        bounds = [(7*pi)/8, pi];
        mult = 2;
    else
        error("invalid diff!");
    end
end

function h = plot_circle(angles, values, g1, g2)
    hold on
    
    rectangle('Position',[-1 -1 2 2],'Curvature',1, 'EdgeColor',[0.5,0.5,0.5]);

    for i=1:length(angles)
        th = angles(i);
        xunit = cos(th);
        yunit = sin(th);
        scatter(xunit, yunit, max(5,values(i)*100), 'k', 'filled');
        if xunit > 0 && yunit <= 0
            t = text(xunit+0.2, yunit-0.1, num2str(values(i), 3));
        elseif xunit >= 0 && yunit > 0
            t= text(xunit+0.2, yunit+0.1, num2str(values(i), 3));
        elseif xunit < 0 && yunit <= 0
            t= text(xunit-0.3, yunit-0.1, num2str(values(i), 3));
        elseif xunit <= 0 && yunit > 0
            t= text(xunit-0.3, yunit+0.1, num2str(values(i), 3));
        end
        t.FontSize = 8;
        c = min(1,values(i)*100);
        c = abs(0.92-c);
        quiver(0, 0, xunit, yunit, 'Color', [c,c,c]);
    end
    scatter(g1(1), g1(2), 100, 'r', 'filled');
    scatter(g2(1), g2(2), 100, 'b', 'filled');
    scatter(0, 0, 50, 'k');
    text(g1(1)+0.1, g1(2),'g1', 'Color', 'r');
    text(g2(1)+0.1, g2(2),'g2', 'Color', 'b');
    
    xlim([-1.5,1.5]);
    ylim([-1.5,1.5]);
    set(gcf,'Position',[100 100 500 500]);
    set(gcf,'color','w');
    hold off
end