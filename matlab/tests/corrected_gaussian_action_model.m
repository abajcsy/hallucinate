clear all
close all
clf

% control bounds.
urange = [-pi, pi];

% goal locations (meters).
g1 = [1; -1];
g2 = [1; 1];

% initial state.
x = [0;0];

% number of discrete controls.
num_ctrls = 61;

% get discretized controls.
us = gen_controls(urange, num_ctrls);

% get control bounds for each control.
ubounds = gen_ubounds(us, urange, num_ctrls);

% standard deviation of gaussian.
mu = 0;
sigma = pi/8;

% truncated gaussian with zero mean.
pd = makedist('Normal','mu',mu,'sigma',sigma);
truncpd = truncate(pd, urange(1), urange(2));

% 50/50 prior on each goal.
priorg1 = 0.5;

% action probabilities
us_probs_g1 = zeros(1,length(us));
us_probs_g2 = zeros(1,length(us));

% posteriors
posteriors = zeros(1,length(us));

% compute optimal action.
uopt_g1 = atan2(g1(2)- x(2), g1(1) - x(1));
uopt_g2 = atan2(g2(2)- x(2), g2(1) - x(1));
    
% get the probability of each action.
for i=1:num_ctrls
    u = us(i);
    
    pu_g1 = compute_prob(u, uopt_g1, us, ubounds, truncpd);
    pu_g2 = compute_prob(u, uopt_g2, us, ubounds, truncpd);
    fprintf(strcat("u: ", num2str(u), ", pu_g1=", num2str(pu_g1), "\n"));
    
    % compute P(g1|u,x) \propto P(u|x,g1)*P(g1)
    posteriorpug1 = (pu_g1 * priorg1)/(pu_g1 * priorg1 + pu_g2 *(1-priorg1));
    
    % plotting info...
    us_probs_g1(i) = pu_g1;
    us_probs_g2(i) = pu_g2;
    posteriors(i) = posteriorpug1;
end

% Tolerance for when computing the argmax_u Posterior(g1 | u). 
tol = 1e-3;

% Plot action probabilities: P(u|x,g1)
figure(1)
plot_circle(us, us_probs_g1, g1, g2, x, posteriors, tol);
title("P(u|x,g1)");

% Plot posterior: P(g1|u,x)
figure(2)
plot_posterior(posteriors, us, g1, x, priorg1, us_probs_g1, tol);

%% Compute probability.
function pu = compute_prob(u, uopt, us, ubounds, truncpd)

    % minimum angular distance between current control (u) and uopt
    diff = abs(angdiff(u, uopt));
    
    % corner case.
    if length(us) == 1
        pu = 1;
        return
    end

    % note because of numerics: sometimes we get controls that are just close to zero but are negative
    zero_tol = -1e-7; 
    % find indicies of all controls in the positive [0, pi] range.
    pos_idxs = find(us >= zero_tol);
    pos_idxs(end+1) = 1; %include -pi since its = pi
    
    % find the control bounds.
    for i=pos_idxs
        bounds = ubounds{i};
        low_bound = bounds(1);
        up_bound = bounds(2);
                
        if low_bound < 0 && diff <= up_bound
            % catch corner case around 0.
            p = cdf(truncpd, [0, up_bound]);
            pu = abs(p(2) - p(1)) * 2;
        elseif up_bound < 0 && diff >= low_bound
            % considering an action that falls into the bounds around -pi
            p = cdf(truncpd, [low_bound, pi]);
            pu = abs(p(2) - p(1)) * 2;
        elseif diff <= up_bound && diff >= low_bound
            % normal integration over bounds.
            p = cdf(truncpd, [low_bound, up_bound]);
            pu = abs(p(2) - p(1));
        end
    end
    
end

%% Gets the control bounds for integration. 
function ubounds = gen_ubounds(us, urange, num_ctrls)
    ubounds = cell(1,num_ctrls);
    incr = (urange(2) - urange(1))/num_ctrls;
    for i=1:num_ctrls
        u = us(i);
        ubounds{i} = [wrapToPi(u - incr/2), wrapToPi(u + incr/2)];
    end
end

%% Generate discrete num_ctrls ranging from [urange(1), urange(2)) 
function us = gen_controls(urange, num_ctrls)
    incr = (urange(2) - urange(1))/num_ctrls;
    us = zeros(1,num_ctrls);
    u = urange(1);
    for i=1:num_ctrls
        us(i) = u;
        u = u + incr;
    end
end

%% Plotting function.
function h = plot_circle(angles, values, g1, g2, x0, posteriors, tol)
    hold on
    
    rectangle('Position',[x0(1)-1 x0(2)-1 2 2],'Curvature',1, 'EdgeColor',[0.5,0.5,0.5]);
    
    % find indicies of max probability value.
    max_prob_u_idx = find(values == max(values));
    
    % find indicies of the max POSTERIOR value in direction of g1.
    max_post_g1 = max(posteriors);
    close_enough_to_max_idx = find(abs(posteriors - max_post_g1) <= tol);
    
    for i=1:length(angles)
        th = angles(i);
        xunit = x0(1) + cos(th);
        yunit = x0(2) + sin(th);
        scatter(xunit, yunit, max(5,values(i)*100), 'k', 'filled');
        %t = text(xunit, yunit, num2str(values(i), 3));
%         if xunit > 0 && yunit <= 0
%             t = text(xunit+0.2, yunit-0.1, num2str(values(i), 3));
%         elseif xunit >= 0 && yunit > 0
%             t= text(xunit+0.2, yunit+0.1, num2str(values(i), 3));
%         elseif xunit < 0 && yunit <= 0
%             t= text(xunit-0.3, yunit-0.1, num2str(values(i), 3));
%         elseif xunit <= 0 && yunit > 0
%             t= text(xunit-0.3, yunit+0.1, num2str(values(i), 3));
%         end
        t.FontSize = 6;
        c = min(1,values(i)*100);
        c = abs(0.92-c);
        if any(max_prob_u_idx == i) && any(close_enough_to_max_idx == i)
            % plot control that is towards goal AND for posterior in cyan
            quiver(x0(1), x0(2), cos(th), sin(th), 'Color', [0,1,1]);
        elseif any(max_prob_u_idx == i)
            % plot optimal control towards goal in red.
            quiver(x0(1), x0(2), cos(th), sin(th), 'Color', [1,0,0]);
        elseif any(close_enough_to_max_idx == i)
            % plot optimal control for posterior in green.
            quiver(x0(1), x0(2), cos(th), sin(th), 'Color', [0,1,0]);
        else
            quiver(x0(1), x0(2), cos(th), sin(th), 'Color', [c,c,c]);
        end
    end
    scatter(g1(1), g1(2), 100, 'r', 'filled');
    scatter(g2(1), g2(2), 100, 'b', 'filled');
    scatter(x0(1), x0(2), 50, 'k');
    text(g1(1)+0.1, g1(2),'g1', 'Color', 'r');
    text(g2(1)+0.1, g2(2),'g2', 'Color', 'b');
    
    xlim([-3,3]);
    ylim([-3,3]);
    set(gcf,'Position',[100 100 600 600]);
    set(gcf,'color','w');
    hold off
end

function plot_posterior(posteriors, us, g1, x, priorg1, us_probs_g1, tol)
    hold on
    plot(us, posteriors, 'mo-')
    plot([-pi, pi], [priorg1, priorg1], 'k--');
    uopt_g1 = atan2(g1(2)- x(2), g1(1) - x(1));
    if uopt_g1 <= 0
        xticks([us(1), uopt_g1, 0, us(end)]);
        xticklabels({num2str(us(1)), num2str(uopt_g1), '0', num2str(us(end))});
    else
        xticks([us(1), 0, uopt_g1, us(end)]);
        xticklabels({num2str(us(1)), '0', num2str(uopt_g1), num2str(us(end))});
    end

    % find indicies of action that is most likely under action model.
    max_prob_u_idx = find(us_probs_g1 == max(us_probs_g1));

    % find indicies of the action that maximally changes POSTERIOR value in direction of g1.
    max_post_g1 = max(posteriors);
    close_enough_to_max_idx = find(abs(posteriors - max_post_g1) <= tol);
    
    % find indicies of the action that moves towards goal AND the
    % posterior.
    both_idx = intersect(max_prob_u_idx, close_enough_to_max_idx);

    for i=1:length(close_enough_to_max_idx)
        % plot line at  optimal actions that shift posterior.
        idx = close_enough_to_max_idx(i);
        plot([us(idx),us(idx)], [0,1], 'g-');
    end
    for i=1:length(max_prob_u_idx)
        % plot line at optimal actions to goal1
        idx = max_prob_u_idx(i);
        plot([us(idx),us(idx)], [0,1], 'r-');
    end
    for i=1:length(both_idx)
        % plot line at  optimal actions that shift posterior AND towards g1
        idx = both_idx(i);
        plot([us(idx),us(idx)], [0,1], 'c-');
    end


    set(gcf,'Position',[100 100 600 600]);
    set(gcf,'color','w');
    grid on
    title("Posterior for G1")
    hold off
end