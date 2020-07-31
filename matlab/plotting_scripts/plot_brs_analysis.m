clear all
clf

%% load data.
load('beta_1_to_0_brs.mat');
goalSetRad = 0;

%% Plotting params.
fontsize = 20;
ticksize = 20;
dt = 0.1;

%% Saving info.
makeVideo = false;

%% Make figure
if makeVideo
    %f = figure('Position', [0, 0, 1000, 1000]);
    figure(1);
    set(gcf,'Color','w');
    hold on

    %% Plot the goal in 2D.
    gtemp2D = createGrid([-8,-8], [8,8], [51, 51]);
    g1Data2D = shapeSphere(gtemp2D, theta, 0.1);
    [~, goalContour] = contour(gtemp2D.xs{1}, gtemp2D.xs{2}, g1Data2D, [0,0.01]);
    goalContour.LineWidth = 6;
    goalContour.EdgeColor = 'r';

    alpha = 0.5;
    extraArgs.sliceDim = 0;
    extraArgs.applyLight = true;
    extraArgs.alpha = alpha;
    color = [7, 71, 107]/255.;
    level = 0;

    % plot initial condition
    initCond = scatter3(x0(1), x0(2), x0(3));
    initCond.MarkerFaceColor = 'k';
    initCond.MarkerEdgeColor = 'k';
    initCond.SizeData = 50;
    grid on;
    box on;
    %box off;

    % Flip the order of the value functions to match the time. This way we
    % start from the beginning of time: V(x,0) and go to V(x,T), just like
    % time.
    valueFunsFlip = flip(valueFuns,4);


    % make a video!
    path = '/home/abajcsy/Documents/Conferences/';
    videoFilename = strcat(path, 'beta_1_to_0_brs.avi');
    vout = VideoWriter(videoFilename,'Motion JPEG AVI');
    vout.Quality = 100;
    vout.FrameRate = 5;
    vout.open;

    for t = 1:length(times)
        h = visSetIm(g, valueFunsFlip(:,:,:,t), color, level, extraArgs);
        h.FaceAlpha = alpha;
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        set(gca, 'ztick', [0, 0.5, 1]);
        set(gca, 'xlim', [-6.8, 6.8]);
        set(gca, 'ylim', [-6.8, 6.8]);
        set(gca, 'zlim', [-0.1,1.1]);
        set(gca, 'fontsize', ticksize);

        titleStr = strcat('$BRS(z, t=-', num2str(t*dt),')$');
        xlabel('$h^x$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
        ylabel('$h^y$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
        zlabel('$b^t(\beta = 1)$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
        title(titleStr, 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')

        %pause(0.01);
        current_frame = getframe(gcf); %gca does just the plot
        writeVideo(vout,current_frame);

        if t ~= length(times)
            delete(h);
        end
    end

    vout.close;
end

% Plot the optimal trajectory.
plotTraj(traj, traj_tau, theta, grid_min, grid_max, goalSetRad, fontsize,ticksize);

%% Plots the state trajectory.
function plotTraj(traj, traj_tau, theta, ...
    grid_min, grid_max, goalSetRad, fontsize, ticksize)
    f2 = figure(2);
    hold on

    % Setup colors.
    startColor = [0,0,0]; %[79., 0., 128.]/255.;
    endColor = [0, 152, 194]/255.; %[255., 143., 255.]/255.;
    r = linspace(startColor(1), endColor(1), length(traj_tau));
    g = linspace(startColor(2), endColor(2), length(traj_tau));
    b = linspace(startColor(3), endColor(3), length(traj_tau));

    % Record state.
    % Plot first point.
    color = [r(1), g(1), b(1)];
    xcurr = traj(1:3, 1);
    plot3(xcurr(1), xcurr(2), xcurr(3), '-o', 'color', color, ...
        'markeredgecolor', color, 'markerfacecolor', color);
    
    % Add first timestamp.
    %txt = strcat('t=', num2str(traj_tau(1)), ', P(\beta=1)=', num2str(xcurr(3)));
    txt = strcat('$b^0(\beta=1)=', num2str(xcurr(3)), '$');
    tp = text(xcurr(1)+0.3, xcurr(2)+0.05, xcurr(3)+0.05, txt);
    tp.Color = color;
    tp.Interpreter = 'Latex';
    tp.FontSize = fontsize-8;
    for t=2:length(traj_tau)
        xprev = traj(1:3, t-1);
        xcurr = traj(1:3, t);
        % Plot point and connection between pts.
        color = [r(t), g(t), b(t)];
        p = plot3([xprev(1), xcurr(1)], [xprev(2), xcurr(2)], [xprev(3), xcurr(3)], '-o', ...
                'Color', color, ...
                'markeredgecolor', color, ...
                'markerfacecolor', color);
        p.LineWidth = 2;
    end
    xcurr = traj(1:3, end);
    %add timestamps
    %txt = strcat('t=', num2str(traj_tau(end)), ', p=', num2str(xcurr(3)));
    txt = strcat('$b^t(\beta=1)=', num2str(xcurr(3)), '$');
    tp = text(xcurr(1)+0.3, xcurr(2)+0.05, xcurr(3)+0.05, txt);
    tp.Color = color;
    tp.Interpreter = 'Latex';
    tp.FontSize = fontsize-8;

    % Plot goals (red is ground truth, grey is other goal).
    rectangle('Position',[theta(1)-goalSetRad ...
                          theta(2)-goalSetRad ...
                          goalSetRad*2 ...
                          goalSetRad*2],...
                          'Curvature',1, ...
                          'FaceColor',[1, 0.67, 0.67],...
                          'EdgeColor',[1, 0.67, 0.67],...
                          'LineWidth',1);
    plot3(theta(1), theta(2), 0.5, '-o', ...
                'Color', 'r', ...
                'markeredgecolor', 'r', ...
                'markerfacecolor', 'r');
    %g1Txt = strcat('g', num2str(1));
    %t1 = text(theta(1), theta(2), 0.55, g1Txt);
    %t1.FontSize = 12;
    %t1.Color = 'r';
    
    grid off;
    box on;
    xlim([grid_min(1), grid_max(1)]);
    ylim([grid_min(2), grid_max(2)]);
    zlim([grid_min(3), grid_max(3)]);
    
    xlabel('$h^x$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    ylabel('$h^y$', 'interpreter', 'latex', 'fontsize', fontsize, 'fontweight','bold')
    set(gca, 'xtick', [traj(1,end), traj(1,1)]);
    set(gca, 'ytick', [traj(2,1), traj(2,end)]);
    set(gca, 'ztick', [0, 0.5, 1]);
    set(gca, 'xlim', [-6.8, 6.8]);
    set(gca, 'ylim', [-6.8, 6.8]);
    set(gca, 'zlim', [-0.1,1.1]);
    %set(gca, 'fontsize', ticksize);
    f2.Position = [0,0,650,650];
    set(gcf,'Color','w');
end