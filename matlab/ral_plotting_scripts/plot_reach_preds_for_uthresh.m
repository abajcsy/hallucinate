% Setup which preds we wanna see.
pg1 = 0.5;
delta = 0.005; %0.005, 0.01, 0.03
saveFig = false;

% Load data.
repo = what('hallucinate');
folder = '/ral_data_uthresh/';
filename = strcat('reach_static_pg1',num2str(pg1),...
    '_gamma0.5_delta',num2str(delta),'_1pred.mat');
load(strcat(repo.path, folder, filename));

totalTime = length(all_preds);

% Known goal locations (in m). 
goals = {[2, 2], [2, -2]}; 

% Grid
grid_min = [-4; -4];  % Lower corner of computation domain
grid_max = [4; 4];     % Upper corner of computation domain

preds = all_preds{1};
xcurr = human_states{1};

% [reachability] Color setup.
reach_start_color = [97, 17, 90]/255.;
reach_end_color = [235, 19, 216]/255.;
r_red = linspace(reach_start_color(1), reach_end_color(1), length(preds(1,1,1,:)));
r_green = linspace(reach_start_color(2), reach_end_color(2), length(preds(1,1,1,:)));
r_blue = linspace(reach_start_color(3), reach_end_color(3), length(preds(1,1,1,:)));
alpha = linspace(0.8,0.6,length(preds(1,1,1,:)));

figure(1)
hold on;

dimsToRemove = [0 0 1];

% Plot the predictions.
predIncr = 2;
small_dt = 0.1667;
vel = 0.6;
for pt=1:predIncr:length(preds(1,1,1,:))
    curr_pred_color = [r_red(pt),  r_green(pt), r_blue(pt)];
   
    % Plot full forward reachable set (FRS) = dt*vel
    DT = (pt-1)*small_dt;
    frs_rad = DT*vel;
    frs_color = [alpha(pt), alpha(pt), alpha(pt)]; %curr_pred_color;
    frs = rectangle('Position',[xcurr(1)-frs_rad xcurr(2)-frs_rad frs_rad*2 frs_rad*2],...
        'Curvature',1, 'EdgeColor', frs_color);
    frs.LineStyle = '--';

    % Plot prediction contour.
    [g2d, data2D] = proj(g, preds(:,:,:,pt), dimsToRemove, 'min');
    [M, c] = contour(g2d.xs{1}, g2d.xs{2}, data2D, [0,0]);
    c.LineWidth = 2;
    c.EdgeColor = curr_pred_color;
end

% Plot state of human.
scatter(xcurr(1), xcurr(2), 60, 'k', 'filled');

% Plot goals.
figure(1);
scatter(goals{1}(1), goals{1}(2), 100, 'r', 'filled');
scatter(goals{2}(1), goals{2}(2), 100, 'r', 'filled');

% Label goal 2.
tg1 = text(goals{1}(1)+0.1, goals{1}(2), '$g_1$', 'Color', 'r', 'Interpreter', 'Latex');
tg1.FontSize = 18;

% Label goal 2.
tg2 = text(goals{2}(1)+0.1, goals{2}(2), '$g_2$', 'Color', 'r', 'Interpreter', 'Latex');
tg2.FontSize = 18;

xlim([-2.4, 2.4]);
ylim([-2.4, 2.4]);
set(gcf,'Position',[100 100 400 400]);
set(gcf,'color','w');
whitebg('w');

box on
set(gca,'xcolor','k','ycolor','k', ...
    'xtick',[], 'xticklabel',[], ...
    'ytick',[], 'yticklabel',[])
set(gcf,'color','w');

if saveFig
    repo = what('hallucinate');
    filename = strcat('preds_pg1', num2str(pg1), '_delta' , ...
        num2str(delta), '.png');
    saveas(gcf, strcat(repo.path, '/ral_imgs_uthresh/', filename));
end