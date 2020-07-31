clear all
%clf
%% Load data.
repo = what('hallucinate');
folder = '/ral_data/';

filename = 'reach_static_pg10.5_gamma0.5_delta0.01_1pred.mat';
reach_start_color = [224, 132, 198]/255.;
reach_end_color = [222, 0, 181]/255.;

goals = {[2;2], [2;-2]};
dt = 0.1667;

load(strcat(repo.path, folder, filename));

%% Plot!
f1 = figure(1);
hold on;

t = 1;

% Plot the predictions.
preds = all_preds{t};
xcurr = human_states{t};

finalTime = length(preds(1,1,1,:));
predIncr = 3;
hcolor = 0.1;
predColorR = linspace(hcolor, reach_end_color(1), finalTime);
predColorG = linspace(hcolor, reach_end_color(2), finalTime);
predColorB = linspace(hcolor, reach_end_color(3), finalTime);

alpha = 0.3;

%% plot initial condition.
pt = 1;
curr_pred_color = [predColorR(pt),  predColorG(pt), predColorB(pt)];
applyLight = true;
visSetIm3D(g, preds(:,:,:,pt), curr_pred_color, 0, alpha, applyLight);

%% Plot action distribution.
% control bounds.
urange = [-pi, pi];

% number of discrete controls.
num_ctrls = 31;

% get discretized controls.
us = gen_controls(urange, num_ctrls);

% standard deviation of gaussian.
mu = 0;
sigma = pi/4;

% truncated gaussian with zero mean.
pd = makedist('Normal','mu',mu,'sigma',sigma);
truncpd = truncate(pd, urange(1), urange(2));

% for state x0: compute optimal action.
uopt_g1_x0 = atan2(goals{1}(2)- xcurr(2), goals{1}(1) - xcurr(1));
uopt_g2_x0 = atan2(goals{2}(2)- xcurr(2), goals{2}(1) - xcurr(1));

prior = [0.5, 0.5];
uThresh = 0.01;

for i=1:num_ctrls
    u = us(i);

    pu_g1 = compute_prob_normalized(u, uopt_g1_x0, us, sigma, urange);
    pu_g2 = compute_prob_normalized(u, uopt_g2_x0, us, sigma, urange);
    
    if abs(u-0.7094) < 0.001 
        q = quiver(xcurr(1), xcurr(2), dt*4*cos(u), dt*4*sin(u), 'Color', 'r');
        q.LineWidth = 2;
        q.MaxHeadSize = 0.5;
        
        xnext = [xcurr(1) + dt*cos(u); xcurr(2) + dt*sin(u)]; 
        posterior = pu_g1*prior(1)/(pu_g1*prior(1) + pu_g2*prior(2));
        phi = atan2(posterior-prior(1), xnext(2)-xcurr(2));
        
        q3 = quiver3(xcurr(1), xcurr(2), prior(1), ...
            dt*2*cos(u), dt*2*sin(u), dt*2*sin(phi), ...
            'Color', 'r');
        q3.LineWidth = 2;
        q3.MaxHeadSize = 1;
        
        next_joint_state = [xcurr(1)+dt*1.9*cos(u), ...
            xcurr(2)+dt*1.9*sin(u), ...
            prior(1)+dt*1.9*sin(phi)];
        scatter3(next_joint_state(1), ...
            next_joint_state(2), ...
            next_joint_state(3), 'r', 'filled');
        %scatter(xnext(1), xnext(2), 'r', 'filled');
        
        txt = strcat('$(', num2str(next_joint_state(1), 2), ...
            ',',num2str(next_joint_state(2), 2), ...
            ',',num2str(next_joint_state(3), 2), ')$');
        tinit = text(next_joint_state(1)+0.1, ...
            next_joint_state(2), ...
            next_joint_state(3), txt);
        tinit.Interpreter = 'Latex';
        tinit.FontSize = 14;
        tinit.Color = 'r';
        
    elseif abs(u+0.7094) < 0.001 
        q = quiver(xcurr(1), xcurr(2), dt*4*cos(u), dt*4*sin(u), 'Color', 'b');
        q.LineWidth = 2;
        q.MaxHeadSize = 0.5;
        
        xnext = [xcurr(1) + dt*cos(u); xcurr(2) + dt*sin(u)]; 
        posterior = pu_g1*prior(1)/(pu_g1*prior(1) + pu_g2*prior(2));
        phi = atan2(posterior-prior(1), xnext(2)-xcurr(2));
        
        q3 = quiver3(xcurr(1), xcurr(2), prior(1), ...
            dt*2*cos(u), dt*2*sin(u), dt*2*sin(phi), ...
            'Color', 'b');
        q3.LineWidth = 2;
        q3.MaxHeadSize = 1;
        
        next_joint_state = [xcurr(1)+dt*1.9*cos(u), ...
            xcurr(2)+dt*1.9*sin(u), ...
            prior(1)+dt*1.9*sin(phi)];
        scatter3(next_joint_state(1), ...
            next_joint_state(2), ...
            next_joint_state(3), 'b', 'filled');
        
        txt = strcat('$(', num2str(next_joint_state(1), 2), ...
            ',',num2str(next_joint_state(2), 2), ...
            ',',num2str(next_joint_state(3), 2), ')$');
        tinit = text(next_joint_state(1)+0.1, ...
            next_joint_state(2), ...
            next_joint_state(3), txt);
        tinit.Interpreter = 'Latex';
        tinit.FontSize = 14;
        tinit.Color = 'b';
    else
    
        if pu_g1*prior(1) + pu_g2*prior(2) >= uThresh
            arrow_color = [0.8,0.8,0.8];
            q = quiver(xcurr(1), xcurr(2), dt*4*cos(u), dt*4*sin(u), 'Color', arrow_color);
            q.LineWidth = 2;
            q.MaxHeadSize = 0.5;
        end
    
    end
    
end

%% Plot 2D initial condition.
dimsToRemove = [0 0 1];
[g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
[M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
c.LineWidth = 3;
c.EdgeColor = curr_pred_color;
Xinner = M(1,2:end);
Yinner = M(2,2:end);
c = fill(Xinner,Yinner, curr_pred_color, ...
    'FaceAlpha', 1);

% label initial cond.
tinit = text(xcurr(1)-1, xcurr(2)+0.2, prior(1), '$(0,0,0.5)$');
tinit.Interpreter = 'Latex';
tinit.FontSize = 14;

%% Plot goals.
scatter3(goals{1}(1), goals{1}(2), 0, 100, 'r', 'filled');
scatter3(goals{2}(1), goals{2}(2), 0, 100, 'r', 'filled');

% Label goal 1.
tg1 = text(goals{1}(1)+0.3, goals{1}(2)-0.1, '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.3, goals{2}(2)+0.1, '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-1, 3]);
ylim([-2.5, 2.5]);
zlim([0, 1]);
xlabel('$h_x$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$h_y$', 'Interpreter', 'Latex', 'FontSize', 20);
zlabel('$b^0(\lambda = g_1)$', 'Interpreter', 'Latex', 'FontSize', 20);
set(gcf,'Position',[100 100 700 600]);
set(gcf,'color','w');
whitebg('w');

title('$\mathcal{V}(0) = \{z : V(0,z) \leq 0\}$', ...
    'Interpreter', 'Latex', 'FontSize', 20);

box off, 
ax = gca
grid(ax, 'on');
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[-1, 0, 1, 2], 'xticklabel',[-1, 0, 1, 2], ...
    'ytick',[-2, -1, 0, 1, 2], 'yticklabel',[-2, -1, 0, 1, 2], ...
    'ztick', [0, 0.5, 1], 'zticklabel', [0,0.5,1]);
set(gcf,'color','w');

%% 3D Visualization
function h = visSetIm3D(g, data, color, level, alpha, applyLight)
% h = visSetIm3D(g, data, color, level, applyLight)
% Visualizes a 3D reachable set

[ mesh_xs, mesh_data ] = gridnd2mesh(g, data);

h = patch(isosurface(mesh_xs{:}, mesh_data, level));
isonormals(mesh_xs{:}, mesh_data, h);
h.FaceColor = color;
h.EdgeColor = 'none';
h.FaceAlpha = alpha;

if applyLight
  lighting phong
  camlight left
  %camlight right
end

view(3)
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
