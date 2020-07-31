clear all
clf
%% Load data.
repo = what('hallucinate');
folder = '/ral_data/';

%filename = 'reach_static_pg10.9_gamma0.5_delta0.01_1pred.mat';
filename = 'reach_static_pg10.5_gamma0.5_delta0.01_1pred.mat';
reach_start_color = [224, 132, 198]/255.;
reach_end_color = [222, 0, 181]/255.;

goals = {[2;2], [2;-2]};

load(strcat(repo.path, folder, filename));

%% Plot!
f1 = figure(2);
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

alpha = 0.2;

for pt=1:predIncr:finalTime
    curr_pred_color = [predColorR(pt),  predColorG(pt), predColorB(pt)];
    
    applyLight = true;
    visSetIm3D(g, preds(:,:,:,pt), curr_pred_color, 0, alpha, applyLight);

    dimsToRemove = [0 0 1];
    [g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
    [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
    c.LineWidth = 3;
    c.EdgeColor = curr_pred_color;
end

scatter3(goals{1}(1), goals{1}(2), 0, 100, 'r', 'filled');
scatter3(goals{2}(1), goals{2}(2), 0, 100, 'r', 'filled');

% Label goal 1.
tg1 = text(goals{1}(1)+0.3, goals{1}(2)-0.1, '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.3, goals{2}(2)+0.1, '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-1.5, 2.5]);
ylim([-2.5, 2.5]);
zlim([0, 1]);
xlabel('$h_x$', 'Interpreter', 'Latex', 'FontSize', 20);
ylabel('$h_y$', 'Interpreter', 'Latex', 'FontSize', 20);
zlabel('$b^{\tau}(\lambda = g_1)$', 'Interpreter', 'Latex', 'FontSize', 20);
set(gcf,'Position',[100 100 700 600]);
set(gcf,'color','w');
whitebg('w');

title('$\mathcal{V}(\tau) = \{z : V(\tau,z) \leq 0\}$', ...
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
