clear all
%close all
clf

% control bounds.
urange = [-pi, pi];

% goal locations (meters).
g1 = [1; -1];
g2 = [1; 1];

% initial state.
x0 = [0; 0];

% number of discrete controls.
num_ctrls = 31;

% get discretized controls.
us = gen_controls(urange, num_ctrls);

% get control bounds for each control.
ubounds = gen_ubounds(us, urange, num_ctrls);

% standard deviation of gaussian.
mu = 0;
sigma = pi/4;

% truncated gaussian with zero mean.
pd = makedist('Normal','mu',mu,'sigma',sigma);
truncpd = truncate(pd, urange(1), urange(2));

% for state x0: compute optimal action.
uopt_g1_x0 = atan2(g1(2)- x0(2), g1(1) - x0(1));
uopt_g2_x0 = atan2(g2(2)- x0(2), g2(1) - x0(1));

% action probabilities
us_probs_g1_x0 = zeros(1,length(us));
us_probs_g2_x0 = zeros(1,length(us));

% get the probability of each action.
for i=1:num_ctrls
    u = us(i);
    
    % Compute probabilities for u and x0 with INTEGRATION scheme.
    %pu_g1 = compute_prob(u, uopt_g1_x0, us, ubounds, truncpd);
    %pu_g2 = compute_prob(u, uopt_g2_x0, us, ubounds, truncpd);
    
    % Compute probabilities for u and x0 with NORMALIZATION scheme.
    pu_g1 = compute_prob_normalized(u, uopt_g1_x0, us, sigma, urange);
    pu_g2 = compute_prob_normalized(u, uopt_g2_x0, us, sigma, urange);

    fprintf(strcat("u: ", num2str(u), ", pu_g1=", num2str(pu_g1), "\n"));
    
    % plotting info...
    us_probs_g1_x0(i) = pu_g1;
    us_probs_g2_x0(i) = pu_g2;
end

% Likley control threshold. 
threshByPercentile = true;
% uThresh = 0.03;
% uThresh = 0.1;
% uThresh = 0.2;
% uThresh = 0.3;
% uThresh = 0.4;
% uThresh = 0.5;
% uThresh = 0.6;
% uThresh = 0.7;
% uThresh = 0.8;
uThresh = 0.9;

% Prior.
prior = [0.5, 0.5];

% If we should save figure.
saveFig = false;

%% Plot state-dependant action probabilities: P(u|x)
f1 = figure(1);
plot_pu(us, us_probs_g1_x0, us_probs_g2_x0, g1, g2, ...
    x0, uThresh, prior, saveFig, threshByPercentile);

%% Plot action probabilities: P(u|x,g1)
figure(2)
plot_pu_given_g(us, us_probs_g1_x0, g1, g2, x0);
title("P(u|x,g1)");

%% Compute probability by integrating on control bounds.
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
    % Need to make sure we grab bounds around zero.
    for i=1:length(us)
        bound = ubounds{i};
        if bound(2) >= 0 && bound(1) <= 0
            pos_idxs(end+1) = i;
        end
    end
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

%% Computes probabilit by querying Gaussian distribution and normalizes. 
function pu = compute_prob_normalized(u, uopt, us, sigma, urange)

    pd = makedist('Normal', 'mu', uopt, 'sigma', sigma);
    truncpd_g = truncate(pd, urange(1), urange(2));
   
    % Get probability of this action under gaussian
    pu_g = pdf(truncpd_g, u);
    
    % Compute normalizer. 
    normalizer_g = 0.0;
    for i=1:length(us)
        ucurr = us(i);
        normalizer_g = normalizer_g + pdf(truncpd_g, ucurr);
    end
    
    pu = pu_g/normalizer_g;
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

%% Return the indices of the controls of a given percentile of likelihood.
function u_idxs = threshold_to_percentile(u_probs, thresh)
    [sorted_u_probs, all_u_idxs] = sort(u_probs, 'descend');

    prob_so_far = 0.0;
    u_idxs = [];
    for i=1:length(sorted_u_probs)
        prob_so_far = prob_so_far + sorted_u_probs(i);
        if prob_so_far < thresh
            u_idxs(end+1) = all_u_idxs(i);
        else
            break;
        end
    end
end

%% Plotting function: P(u | x) = P(u | x, g1)P(g1) + P(u | x, g2)P(g2)
function h = plot_pu(angles, pug1, pug2, g1, g2, x0, uThresh, prior, saveFig, threshByPercentile)
    hold on
    
    %rectangle('Position',[x0(1)-1 x0(2)-1 2 2],...
    %    'Curvature',1, 'EdgeColor', [0.8,0.8,0.8]);
    
    offset = 7;
    
    % Setup colormap and plotting.
    cmap = colormap(flipud(gray));
    cmap_smaller = cmap(offset:end-offset, :);
    
    pug1_thresholded = pug1;
    pug2_thresholded = pug2;
    
    if threshByPercentile
        pug1_thresholded_uidxs = threshold_to_percentile(pug1, 1 - uThresh);
        pug2_thresholded_uidxs = threshold_to_percentile(pug2, 1 - uThresh);
    
        pug1_thresholded(setdiff(1:length(pug1), pug1_thresholded_uidxs)) = 0;    
        pug2_thresholded(setdiff(1:length(pug2), pug2_thresholded_uidxs)) = 0;
    end
    
    all_probs = pug1_thresholded*prior(1) + pug2_thresholded*prior(2);
    max_prob = max(all_probs);
    
    norm = 1/max_prob;
    normalized_cmap = cmap/norm;
    normalized_cmap = normalized_cmap(offset:end-offset, :);
    
    colormap(cmap_smaller);
    cbar = colorbar;
    caxis([0, max_prob]);
    %caxis([0, 0.2]);
    
    for i=1:length(angles)
        th = angles(i);
        xunit = x0(1) + cos(th);
        yunit = x0(2) + sin(th);
        
        % Compute: P(u | x) = P(u | x, g1)P(g1) + P(u | x, g2)P(g2)
%         prob = pug1(i)*prior(1) + pug2(i)*prior(2);
        prob = pug1_thresholded(i)*prior(1) + pug2_thresholded(i)*prior(2);
        
        if (threshByPercentile && prob > 0) || ...
           (~threshByPercentile && prob >= uThresh)
            if xunit > 0 && yunit <= 0
                t = text(xunit+0.2, yunit-0.1, num2str(prob, 2));
            elseif xunit >= 0 && yunit > 0
                t= text(xunit+0.2, yunit+0.1, num2str(prob, 2));
            elseif xunit < 0 && yunit <= 0
                t= text(xunit-0.4, yunit-0.1, num2str(prob, 2));
            elseif xunit <= 0 && yunit > 0
                t= text(xunit-0.4, yunit+0.1, num2str(prob, 2));
            end
            c = min(0.92,prob*10);
            c = abs(0.92-c);
            t.FontSize = 10;
            t.Interpreter = 'Latex';
            t.Color = [c,c,c];
            
            inverted_prob = (max_prob-prob);
            [minDistance, indexOfMin] = min(abs(normalized_cmap(:,1) - inverted_prob));
            
            curr_color = cmap_smaller(indexOfMin, :);
            q = quiver(x0(1), x0(2), cos(th), sin(th), 'Color', curr_color);
            q.LineWidth = 2;
            q.MaxHeadSize = 0.5;
            
            s = scatter(xunit, yunit, max(3,prob*500));
            s.MarkerEdgeColor = curr_color;
            %s.MarkerFaceColor = curr_color;
        end
    end
    
    % Plot goals and human initial state.
    scatter(g1(1), g1(2), 100, 'r', 'filled');
    scatter(g2(1), g2(2), 100, 'b', 'filled');
    scatter(x0(1), x0(2), 90, 'k', 'filled');
    
    % Label goal 2.
    text_pg1 = strcat("$P^t(g_1) = ", num2str(prior(1)), "$");
    tpg1 = text(g1(1)-0.4, g1(2)-0.3, text_pg1, 'Color', 'r', 'Interpreter', 'Latex');
    tg1 = text(g1(1)+0.1, g1(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
    tg1.FontSize = 15;
    tpg1.FontSize = 14;
    
    % Label goal 2.
    text_pg2 = strcat("$P^t(g_2) = ", num2str(prior(2)), "$");
    tpg2 = text(g2(1)-0.4, g2(2)+0.3, text_pg2, 'Color', 'b', 'Interpreter', 'Latex');
    tg2 = text(g2(1)+0.1, g2(2), '$g_2$', 'Color', 'b', 'Interpreter', 'Latex');
    tg2.FontSize = 15;
    tpg2.FontSize = 14;
    
    % Setup overall figure formatting.
    xlim([-1.5,1.8]);
    ylim([-1.8,1.8]);
    set(gcf,'Position',[100 100 500 500]);
    whitebg('w');
    box on
    set(gca,'xcolor','k','ycolor','k', ...
        'xtick',[], 'xticklabel',[], ...
        'ytick',[], 'yticklabel',[])
    set(gcf,'color','w');
    
    % Plot title.
    tstr = strcat("$\mathcal{U}(z^t) = \{ u : P^t(u^t_H \mid z^t) \geq ",...
    num2str(uThresh), "\}$");
    t = title(tstr, 'Interpreter', 'Latex');
    t.FontSize = 18;
    hold off
    
    if saveFig
        repo = what('hallucinate');
        filename = strcat('pg1', num2str(prior(1)), '_uthr' , ...
            num2str(uThresh), '.png');
        saveas(gcf, strcat(repo.path, '/ral_imgs/', filename));
    end
end

%% Plotting function: P(u | x, g)
function h = plot_pu_given_g(angles, values, g1, g2, x0)
    hold on
    
    rectangle('Position',[x0(1)-1 x0(2)-1 2 2],'Curvature',1, 'EdgeColor',[0.5,0.5,0.5]);
    
    for i=1:length(angles)
        th = angles(i);
        xunit = x0(1) + cos(th);
        yunit = x0(2) + sin(th);
        scatter(xunit, yunit, max(5,values(i)*100), 'k', 'filled');
        %t = text(xunit, yunit, num2str(values(i), 3));
        if xunit > 0 && yunit <= 0
            t = text(xunit+0.2, yunit-0.1, num2str(values(i), 3));
        elseif xunit >= 0 && yunit > 0
            t= text(xunit+0.2, yunit+0.1, num2str(values(i), 3));
        elseif xunit < 0 && yunit <= 0
            t= text(xunit-0.4, yunit-0.1, num2str(values(i), 3));
        elseif xunit <= 0 && yunit > 0
            t= text(xunit-0.4, yunit+0.1, num2str(values(i), 3));
        end
        t.FontSize = 8;
        c = min(1,values(i)*10);
        c = abs(0.92-c);
        q = quiver(x0(1), x0(2), cos(th), sin(th), 'Color', [c,c,c]);
        q.LineWidth = max(3,values(i)*10);
    end
    scatter(g1(1), g1(2), 100, 'r', 'filled');
    scatter(g2(1), g2(2), 100, 'b', 'filled');
    scatter(x0(1), x0(2), 50, 'k');
    text(g1(1)+0.1, g1(2),'g1', 'Color', 'r');
    text(g2(1)+0.1, g2(2),'g2', 'Color', 'b');
    
    xlim([-3,3]);
    ylim([-3,3]);
    set(gcf,'Position',[100 100 800 800]);
    set(gcf,'color','w');
    whitebg('w');
    hold off
end
