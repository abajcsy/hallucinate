clear all
close all

%% Load data.
repo = what('hallucinate');
filename = 'sigma_init_cond.mat';
load(strcat(repo.path, '/ral_data_revise/human_traj_data/', filename));

%% (Bayesian) Predictor params.

% Known human goal locations. 
goals = {[2; 2], [2; -2]}; 


%% Pick which sigma and init condition to debug.
sigma = human_sigmas(5);
x0 = human_init_conds{1};

%% Run metrics!

traj_map = all_trajs(num2str(sigma));
human_traj = traj_map(num2str(x0));

hold on
for t=1:simT
    xcurr = human_traj{t};
    val = (simT - t)/simT;
    color = [val, val, val];
    % plot true goal.
    scatter(goals{1}(1), goals{1}(2), 99, 'r', 'filled');
    
    s = scatter(xcurr(1), xcurr(2), 40, color, 'filled');
    s.MarkerEdgeColor = 'k';
    s.LineWidth = 1;
    xlim([-4, 4]);
    ylim([-4, 4]);
    grid on
    set(gcf,'Position',[100 100 700 700]);
    set(gcf,'color','w');
end
