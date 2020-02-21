clf
clear all
close all

%load('northwestern_frs_5050_vis.mat');
%load('northwestern_frs_9010_vis.mat');
load('northwestern_frs_9010_01thresh_vis.mat');
%load('northwestern_frs_9010_03thresh_vis.mat');

% params = frsTwoGoalPrecompute();
% 
% predictor = HumanPredictor(params);
% predictor.updatePredictions();
% 
% % Planning step.
% [data, times] = predictor.getPredictions();

% Plotting params.
fontsize = 60;
ticksize = 20;

%% Show the value function propagation
path = '/home/abajcsy/Documents/Conferences/ICRA2020_HallucinatingData/';
videoFilename = strcat(path, 'frs_joint_state_9010.avi');
vout = VideoWriter(videoFilename,'Motion JPEG AVI');
vout.Quality = 100;
vout.FrameRate = 5;
vout.open;

f = figure('Position', [0, 0, 1000, 1000]);
set(gcf,'Color','w');
hold on
data_size = size(data);
greyColor = [184, 184, 184]/255.;
magentaColor = [247, 0, 157]/255.;
colorR = linspace(greyColor(1), magentaColor(1), data_size(4));
colorG = linspace(greyColor(2), magentaColor(2), data_size(4));
colorB = linspace(greyColor(3), magentaColor(3), data_size(4));
alpha = linspace(0.5, 0.5, data_size(4));

% Plot the goal 1 in 2D.
gtemp2D = createGrid([-4,-4], [4,4], [81, 81]);
g1Data2D = shapeSphere(gtemp2D, [2;2], 0.1);
[~, c2] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, g1Data2D, [0,0.01]);
c2.LineWidth = 2.5;
c2.EdgeColor = 'r';

% Plot the goal 1 in 2D.
g2Data2D = shapeSphere(gtemp2D, [2;-2], 0.1);
[~, c3] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, g2Data2D, [0,0.01]);
c3.LineWidth = 2.5;
c3.EdgeColor = 'r';

for ii = 1:data_size(4)
    vx = data(:,:,:,ii);  

    % Plot the 3D set. 
    [ mesh_xs, mesh_data ] = gridnd2mesh(predictor.grid, vx);

    h = patch(isosurface(mesh_xs{:}, mesh_data, 0));
    isonormals(mesh_xs{:}, mesh_data, h);
    h.FaceColor = [colorR(ii), colorG(ii), colorB(ii)];
    h.EdgeColor = 'none';
    h.FaceAlpha = alpha(ii);
    
    % Plot the 2D projection.
    [grid2D, data2DFRS] = proj(predictor.grid, vx, [0 0 1], 'min');
    [~, s] = contour(grid2D.xs{1}, grid2D.xs{2}, data2DFRS, [0, 0], ...
        'linecolor', [colorR(ii), colorG(ii), colorB(ii)], 'linewidth', 3);

    lighting phong
    c = camlight;
    c.Position = [100 100 100];
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'ztick', [0, 0.5, 1]);
    set(gca, 'xlim', [-4, 4]);
    set(gca, 'ylim', [-4, 4]);
    set(gca, 'zlim', [-0.1,1.1]);
    set(gca, 'fontsize', ticksize);
    titleStr = strcat('$FRS(z, t=', num2str(ii*params.dt),')$');
    xlabel('$h^x$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    ylabel('$h^y$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    zlabel('$P(g_1)$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    title(titleStr, 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    %zlabel(strcat('$V(', num2str(0.02*(ii-1), '%4.2f'), ', x)$'), 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    view(-58, 16)

    current_frame = getframe(gcf); %gca does just the plot
    writeVideo(vout,current_frame);
    delete(h);
    delete(c);
    delete(s);
end
vout.close;