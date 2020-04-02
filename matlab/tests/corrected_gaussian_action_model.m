clear all
clf

% control bounds.
urange = [-pi, pi];

% goal locations (meters).
g1 = [1; -1];
g2 = [1; 1];

% initial state.
x = [0; 0];

% number of discrete controls.
num_ctrls = 12;

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

% action probabilities
us_probs_g1 = zeros(1,length(us));
us_probs_g2 = zeros(1,length(us));

% get the probability of each action.
for i=1:num_ctrls
    u = us(i);
    
    % comput optimal action.
    uopt_g1 = atan2(g1(2)- x(2), g1(1) - x(1));
    uopt_g2 = atan2(g2(2)- x(2), g2(1) - x(1));
    
    pu_g1 = compute_prob(u, uopt_g1, us, ubounds, truncpd);
    pu_g2 = compute_prob(u, uopt_g2, us, ubounds, truncpd);
    fprintf(strcat("u: ", num2str(u), ", pu_g1=", num2str(pu_g1), "\n"));
    
    % plotting info...
    us_probs_g1(i) = pu_g1;
    us_probs_g2(i) = pu_g2;
end

% Plot action probabilities: P(u|x,g1)
figure(1)
plot_circle(us, us_probs_g1, g1, g2);
title("P(u|x,g1)");

% Plot action probabilities: P(u|x,g2)
%figure(2)
%plot_circle(us, us_probs_g2, g1, g2);
%title("P(u|x,g2)");


%% Compute probability.
function pu = compute_prob(u, uopt, us, ubounds, truncpd)
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
    
    % find the control bounds.
    for i=pos_idxs
        bounds = ubounds{i};
        low_bound = bounds(1);
        up_bound = bounds(2);
        
        if (diff >= low_bound && diff < up_bound)
            if low_bound < 0
                % catch corner case around 0.
                p = cdf(truncpd, [0, up_bound]);
                pu = abs(p(2) - p(1)) * 2;
            else
                % normal integration over bounds.
                p = cdf(truncpd, [low_bound, up_bound]);
                pu = abs(p(2) - p(1));
            end
            break;
        elseif (diff >= low_bound && diff <= pi)
            % considering an action that falls into the bounds around -pi
            p = cdf(truncpd, [up_bound, pi]);
            pu = abs(p(2) - p(1)) * 2;
        end
    end
    
end

%% Gets the control bounds for integration. 
function ubounds = gen_ubounds(us, urange, num_ctrls)
    ubounds = cell(1,num_ctrls);
    incr = (urange(2) - urange(1))/num_ctrls;
    for i=1:num_ctrls
        u = us(i);
        
        if u == -pi || u == pi
            % catch corner case around -pi/pi.
            ubounds{i} = [-wrapToPi(u - incr/2), -wrapToPi(u + incr/2)];
        else
            ubounds{i} = [wrapToPi(u - incr/2), wrapToPi(u + incr/2)]; 
        end
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