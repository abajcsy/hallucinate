clear all
close all
clf

% control bounds.
urange = [-pi, pi];

% goal locations (meters).
g1 = [1; 0];

% initial state.
x = [-4; 0];

v = 0.2;
deltaT = 1;

% number of discrete controls.
num_ctrls = 51;

% get discretized controls.
us = gen_controls(urange, num_ctrls);

% 50/50 prior on each goal.
prior_0 = 0.9;
prior_0_1 = 0;
prior_1 = 0.1;

% action probabilities
us_probs_0 = zeros(1,length(us));
us_probs_0_1 = zeros(1,length(us));
us_probs_1 = zeros(1,length(us));
us_probs = zeros(1,length(us));

% posteriors
posteriors_0 = zeros(1,length(us));
posteriors_0_1 = zeros(1,length(us));
posteriors_1 = zeros(1,length(us));

% get the probability of each action.
for i=1:num_ctrls
    u = us(i);
    
    pu_0 = pu_beta(u,x,0,us,v,deltaT,g1);
    pu_0_1 = pu_beta(u,x,0.1,us,v,deltaT,g1);
    pu_1 = pu_beta(u,x,1,us,v,deltaT,g1);
    pu = pu_0 .* prior_0 + pu_0_1 .* prior_0_1 + pu_1 .* prior_1;
    
    % compute P(g1|u,x) \propto P(u|x,g1)*P(g1)
    posterior_0 = (pu_0 * prior_0)/((pu_0 * prior_0) + (pu_0_1 * prior_0_1) + (pu_1 * prior_1));
    posterior_0_1 = (pu_0_1 * prior_0_1)/((pu_0 * prior_0) + (pu_0_1 * prior_0_1) + (pu_1 * prior_1));
    
    % plotting info...
    us_probs_0(i) = pu_0;
    us_probs_0_1(i) = pu_0_1;
    us_probs_1(i) = pu_1;
    us_probs(i) = pu;
    posteriors_0(i) = posterior_0;
    posteriors_0_1(i) = posterior_0_1;
end

figure(1)
plot_prob(us_probs_0,us,"P(u|x,b=0)");

figure(2)
plot_prob(us_probs_0_1,us,"P(u|x,b=0.1)");

figure(3)
plot_prob(us_probs_1,us,"P(u|x,b=1)");

figure(4)
plot_prob(us_probs,us,"P(u|x)");

% % Tolerance for when computing the argmax_u Posterior(g1 | u). 
% tol = 1e-3;
% 
% % Plot action probabilities: P(u|x,g1)
% figure(1)
% plot_circle(us, us_probs_g1, g1, g2, x, posteriors, tol);
% title("P(u|x,g1)");
% 
% % Plot posterior: P(g1|u,x)
% figure(2)
% plot_posterior(posteriors, us, g1, x, priorg1, us_probs_g1, tol);
% 
% % Plot "continuous-time" posterior: 
% %   (1/dt)*(P(g1|u,x) - P(g1))
% figure(3)
% cont_time_posteriors = (1/deltaT)*(posteriors - priorg1);
% plot_cont_time_posterior(cont_time_posteriors, us, g1, x, priorg1, us_probs_g1, tol);

%% Compute probability (boltzmann)
function pu = pu_beta(u,x,beta,us,v,deltaT,goal)
    us_shape = size(us);
    % Calculta sum of 
    sumControls = 0.0;
    for i = 1:us_shape(2)

        % Get discrete control.
        u_i = us(i);

        % Compute the Q-value of each state and control.
        qval = qFunction(x, u_i, goal, v, deltaT);

        % Calculate value in summation: 
        %   e^{-||(x_t + \Deltat t f(x_t,u_t)) - \theta||_2}
        u_i_val = exp(beta .* qval);

        % Add to running value of summation
        sumControls = sumControls + u_i_val; 
    end
    u_val = exp(beta .* qFunction(x, u, goal, v, deltaT));
    pu = u_val ./ sumControls;
end

function qval = qFunction(x, u, theta, v, deltaT)
    % Find next x by forward euler
    x1 = x(1) + deltaT * v * cos(u);
    x2 = x(2) + deltaT * v * sin(u);

    % Evaluate distance of next x to goal theta under L2 norm
    qval = -1*((x1 - theta(1))^2 + (x2 - theta(2))^2)^(1.);
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

function plot_prob(pu,us,name)
    hold on
    plot(us, pu, 'mo-');
    title(name)
end

%% Plotting function.
function h = plot_circle(angles, values, g1, x0, posteriors, tol)
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
        %scatter(xunit, yunit, max(5,values(i)*100), 'k', 'filled');
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
        c = min(1,values(i)*10);
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
    scatter(x0(1), x0(2), 50, 'k');
    text(g1(1)+0.1, g1(2),'g1', 'Color', 'r');
    
    xlim([-3,3]);
    ylim([-3,3]);
    set(gcf,'Position',[100 100 500 500]);
    set(gcf,'color','w');
    hold off
end

%% Plot discrete posterior update.
function plot_posterior(posteriors, us, g1, x, priorg1, us_probs_g1, tol)
    hold on
    plot(us, posteriors, 'mo-')
    plot([-pi, pi], [priorg1, priorg1], 'k--');
    uopt_g1 = atan2(g1(2)- x(2), g1(1) - x(1));
    if uopt_g1 < 0
        xticks([us(1), uopt_g1, 0, us(end)]);
        xticklabels({num2str(us(1)), num2str(uopt_g1), '0', num2str(us(end))});
    elseif uopt_g1 == 0
        xticks([us(1), uopt_g1, us(end)]);
        xticklabels({num2str(us(1)), num2str(uopt_g1), num2str(us(end))});
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


    set(gcf,'Position',[100 100 500 500]);
    set(gcf,'color','w');
    grid on
    title("Posterior for G1")
    hold off
end

%% Plot continous-time posterior update.
function plot_cont_time_posterior(posteriors, us, g1, x, priorg1, us_probs_g1, tol)
    hold on
    plot(us, posteriors, 'mo-');
    uopt_g1 = atan2(g1(2)- x(2), g1(1) - x(1));
    if uopt_g1 < 0
        xticks([us(1), uopt_g1, 0, us(end)]);
        xticklabels({num2str(us(1)), num2str(uopt_g1), '0', num2str(us(end))});
    elseif uopt_g1 == 0
        xticks([us(1), uopt_g1, us(end)]);
        xticklabels({num2str(us(1)), num2str(uopt_g1), num2str(us(end))});
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
        plot([us(idx),us(idx)], [-0.5,0.5], 'g-');
    end
    for i=1:length(max_prob_u_idx)
        % plot line at optimal actions to goal1
        idx = max_prob_u_idx(i);
        plot([us(idx),us(idx)], [-0.5,0.5], 'r-');
    end
    for i=1:length(both_idx)
        % plot line at  optimal actions that shift posterior AND towards g1
        idx = both_idx(i);
        plot([us(idx),us(idx)], [-0.5,0.5], 'c-');
    end

    set(gcf,'Position',[100 100 500 500]);
    set(gcf,'color','w');
    xlabel("$u \in [-\pi, \pi]$", 'Interpreter','latex');
    ylabel("$\dot{P}_t(g_1)$", 'Interpreter','latex');
    xlim([-pi, pi]);
    ylim([-0.5,0.5]);
    grid on
    title("(Continuous-time) Posterior for G1")
    hold off
end