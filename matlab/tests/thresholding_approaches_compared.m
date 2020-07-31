close all
clear all

% save directory.
folder = ...
    '/home/abajcsy/Documents/Conferences/RAL2020_RobustControlForHumanMotion/threshold_comparisons_v3/';
saveFig = true;

% initial state.  (in grid cell location)
x0 = [4,4];

% right, upright, up, upleft, left, leftdown, down, rightdown
us = {1,2,3,4,5,6,7,8};

% goals: g1 and g2 (in grid cell location)
goals = {[1; 7], [7; 7]};

% prior over goals
b0 = [0.9, 0.1];

% threshold for clipping controls.
delta_action = 0.08;

delta_belief = delta_action*length(goals);

% debugging! print the action probabilities for each goal.
for gIdx=1:length(goals)
    fprintf('------------ g = g%d ----------\n', gIdx);
    for j=1:length(us)
        u = us{j};
        pu = pugiveng(x0, u, gIdx, us, goals);
        fprintf('    u=%d, P(u | x, g) = %f\n', u, pu);
    end
    %fh = plot_ctrls(x0, us, b0, gIdx, goals, 10+gIdx);
end

% compute method via action thresholding.
%preds_action = pred_action(x0, us, goals, b0, delta_action);

% compute method via action thresholding.
preds_action_and_b0 = pred_action_and_b0(x0, us, goals, b0, delta_action);

% compute method via belief tracking thresholding.
preds_belief = pred_belief(x0, us, goals, b0, delta_belief);

% % plotting.
% f1 = figure(1);
% for kn=1:3
%     subplot(1,3,kn);
%     str = strcat('P(x',num2str(kn-1), ' | x0)');
%     h = heatmap(preds_action{kn});
%     h.Colormap = flipud(hot);
%     h.Title = str;
%     h.ColorLimits = [0 0.3];
%     h.CellLabelFormat = '%.2f';
%     h.GridVisible = 'off';
%     if kn == 1 || kn == 2
%         h.ColorbarVisible = 'off';
%     end
% end
% sgtitle(strcat('Action Method: delta_a = ', num2str(delta_action), ', and b0(g1)=', num2str(b0(1))));
% f1.Position = [10 10 1200 400];

% plotting.
f1 = figure(1);
for kn=1:3
    subplot(1,3,kn);
    str = strcat('P(x',num2str(kn-1), ' | x0)');
    h = heatmap(preds_action_and_b0{kn});
    h.Colormap = flipud(hot);
    h.Title = str;
    h.ColorLimits = [0 0.3];
    h.CellLabelFormat = '%.2f';
    h.GridVisible = 'off';
    if kn == 1 || kn == 2
        h.ColorbarVisible = 'off';
    end
end
sgtitle(strcat('Action and Belief Method: delta_a = ', ...
    num2str(delta_action), ', and b0(g1)=', num2str(b0(1))));
f1.Position = [10 10 1200 400];

% plotting.
f2 = figure(2);
for kb=1:3
    subplot(1,3,kb);
    str = strcat('P(x', num2str(kb-1),' | x0)');
    h = heatmap(preds_belief{kb});
    h.Colormap = flipud(hot);
    h.Title = str;
    h.ColorLimits = [0 0.3];
    h.CellLabelFormat = '%.2f';
    h.GridVisible = 'off';
    if kb == 1 || kb == 2
        h.ColorbarVisible = 'off';
    end
end
sgtitle(strcat('Belief Method: delta_b = ', ...
    num2str(delta_belief), ', and b0(g1)=', num2str(b0(1))));
f2.Position = [10 420 1200 400];

% save!
if saveFig
%     filename_action = strcat('action_delta_a = ', ...
%        num2str(delta_action), ', b0(g1)=', num2str(b0(1)));
%     saveas(f1,strcat(folder, filename_action, '.png'));

    filename_action = strcat('action_and_b0_delta_a = ', ...
        num2str(delta_action), ', b0(g1)=', num2str(b0(1)));
    saveas(f1,strcat(folder, filename_action, '.png'));
    
    filename_belief = strcat('belief_delta_b = ', ....
        num2str(delta_belief), ', b0(g1)=', num2str(b0(1)));
    saveas(f2,strcat(folder, filename_belief, '.png'));
end

%% Compute preds via belief thresholding.
function preds_belief = pred_belief(x0, us, goals, b0, delta)
   % size of X and Y dim of state space.
    N = 7;
    
    % Store t=0, t=1, t=2.
    preds_belief = cell(3,1);
    
    %% t=0: store initial state.
    px0 = zeros(N,N);
    px0(x0(1), x0(2)) = 1;
    preds_belief{1} = px0;
    
    %% t=1: store next state. 
    px_t1 = zeros(N,N);

    % Threshold by: \sum_g P(u | x, g)*belief(g) >= delta
    likely_ctrls = belief_thresh_ctrls(x0, us, goals, delta, b0);
 
    for i=1:length(likely_ctrls)
        u = likely_ctrls{i};
        [xnext, ~] = dynamics(x0, u);
        px_t1(xnext(1), xnext(2)) = 1/length(likely_ctrls);
    end
    preds_belief{2} = px_t1;
    
    %% t=2: store next state.
    px_t2 = zeros(N,N);
    
    % Threshold by: \sum_g P(u | x, g)*belief(g) >= delta
    likely_ctrls_t0 = belief_thresh_ctrls(x0, us, goals, delta, b0);
 
    for i0=1:length(likely_ctrls_t0)
        u0 = likely_ctrls_t0{i0};
        x1Idxs = find(px_t1 > 0);
        for iIdx=1:length(x1Idxs)
            [x1_yidx, x1_xidx] = ind2sub(size(px_t1), x1Idxs(iIdx));
            x1 = [x1_yidx, x1_xidx];
            % update belief.
            b1 = next_belief(b0, x0, u0, us, goals);
            % get the new thresholded controls.
            likely_ctrls_t1 = belief_thresh_ctrls(x1, us, goals, delta, b1);
            for i1=1:length(likely_ctrls_t1)
                u1 = likely_ctrls_t1{i1};
                [xnext, ~] = dynamics(x1, u1);
                px_t2(xnext(1), xnext(2)) = 1/length(likely_ctrls_t1) * 1/length(likely_ctrls_t0);
            end
        end
    end
    preds_belief{3} = px_t2;
end

%% Compute preds via action thresholding.
function preds_action = pred_action(x0, us, goals, b0, delta)
    % size of X and Y dim of state space.
    N = 7;
    
    % Store t=0, t=1, t=2.
    preds_action = cell(3,1);
    
    %% t=0: store initial state.
    px0 = zeros(N,N);
    px0(x0(1), x0(2)) = 1;
    preds_action{1} = px0;
    
    %% t=1: store next state. 
    px_t1_goals = cell(2,1);  
    % Compute: P(x1 | x0, g)
    for gIdx=1:length(goals)
        px_t1_g = zeros(N,N);
        
        % Threshold by: P(u | x, g) >= delta
        likely_ctrls = action_thresh_ctrls(x0, us, gIdx, goals, delta);
        
        for i=1:length(likely_ctrls)
            u = likely_ctrls{i};
            [xnext, ~] = dynamics(x0,u);
            % Store {x1, p(x1 | x0, u0, g)}
            px_t1_g(xnext(1), xnext(2)) = 1/length(likely_ctrls);
        end
        
        % put the state probabilities in the list.
        px_t1_goals{gIdx} = px_t1_g;
    end
    % combine to get P(x1 | x0) = \sum_g P(x1 | x0, g) b0(g)
    px_t1 = zeros(N,N);
    for gIdx=1:length(goals)
        px_t1 = px_t1 + px_t1_goals{gIdx} .* b0(gIdx);
    end
    preds_action{2} = px_t1;
    
    %% t=2: store next state. 
    px_t2_goals = cell(2,1); 
    % Compute: P(x2 | x0, g)
    for gIdx=1:length(goals)
        px_t2_g = zeros(N,N);
        for yidx=1:N
            for xidx=1:N
                x = [yidx; xidx];
                
                % optimization!
                if px_t1_goals{gIdx}(yidx, xidx) > 0 
                    % Threshold by: P(u | x, g) >= delta
                    likely_ctrls = action_thresh_ctrls(x, us, gIdx, goals, delta);
                    
                    for i=1:length(likely_ctrls)
                        u = likely_ctrls{i};
                        [xnext, ~] = dynamics(x,u);
                        px_t2_g(xnext(1), xnext(2)) = ...
                            1/length(likely_ctrls) * px_t1_goals{gIdx}(yidx, xidx);
                    end
                end  
            end
        end
        % put the state probabilities in the list.
        px_t2_goals{gIdx} = px_t2_g;
    end 
    % combine to get P(x2 | x0) = \sum_g P(x2 | x0, g) b0(g)
    px_t2 = zeros(N,N);
    for gIdx=1:length(goals)
        px_t2 = px_t2 + px_t2_goals{gIdx} .* b0(gIdx);
    end
    preds_action{3} = px_t2;
end

%% Compute preds via action and b0 thresholding.
function preds_action = pred_action_and_b0(x0, us, goals, b0, delta)
    % size of X and Y dim of state space.
    N = 7;
    
    % Store t=0, t=1, t=2.
    preds_action = cell(3,1);
    
    %% t=0: store initial state.
    px0 = zeros(N,N);
    px0(x0(1), x0(2)) = 1;
    preds_action{1} = px0;
    
    %% t=1: store next state. 
    px_t1_goals = cell(2,1);  
    % Compute: P(x1 | x0, g)
    for gIdx=1:length(goals)
        px_t1_g = zeros(N,N);
        
        % Threshold by: P(u | x, g)*b0(g) >= delta
        likely_ctrls = action_and_b0_thresh_ctrls(x0, us, gIdx, goals, delta, b0);
        
        for i=1:length(likely_ctrls)
            u = likely_ctrls{i};
            [xnext, ~] = dynamics(x0,u);
            % Store {x1, p(x1 | x0, u0, g)}
            px_t1_g(xnext(1), xnext(2)) = 1/length(likely_ctrls);
        end
        
        % put the state probabilities in the list.
        px_t1_goals{gIdx} = px_t1_g;
    end
    % combine to get P(x1 | x0) = \sum_g P(x1 | x0, g) b0(g)
    px_t1 = zeros(N,N);
    for gIdx=1:length(goals)
        px_t1 = px_t1 + px_t1_goals{gIdx} .* b0(gIdx);
    end
    preds_action{2} = px_t1;
    
    %% t=2: store next state. 
    px_t2_goals = cell(2,1); 
    % Compute: P(x2 | x0, g)
    for gIdx=1:length(goals)
        px_t2_g = zeros(N,N);
        for yidx=1:N
            for xidx=1:N
                x = [yidx; xidx];
                
                % optimization!
                if px_t1_goals{gIdx}(yidx, xidx) > 0 
                    % Threshold by: P(u | x, g)*b0(g) >= delta
                    likely_ctrls = action_and_b0_thresh_ctrls(x, us, gIdx, goals, delta, b0);
                    
                    for i=1:length(likely_ctrls)
                        u = likely_ctrls{i};
                        [xnext, ~] = dynamics(x,u);
                        px_t2_g(xnext(1), xnext(2)) = ...
                            1/length(likely_ctrls) * px_t1_goals{gIdx}(yidx, xidx);
                    end
                end  
            end
        end
        % put the state probabilities in the list.
        px_t2_goals{gIdx} = px_t2_g;
    end 
    % combine to get P(x2 | x0) = \sum_g P(x2 | x0, g) b0(g)
    px_t2 = zeros(N,N);
    for gIdx=1:length(goals)
        px_t2 = px_t2 + px_t2_goals{gIdx} .* b0(gIdx);
    end
    preds_action{3} = px_t2;
end

%% Action thresholding. Returns set of likely controls under model.
function likely_ctrls = action_thresh_ctrls(x, us, gIdx, goals, delta)
    likely_ctrls = {};
    for i=1:length(us)
        u = us{i};
        % Threshold by: P(u | x, g) >= delta
        pu = pugiveng(x, u, gIdx, us, goals);
        if pu >= delta
            likely_ctrls{end+1} = u;
        end
    end
end

%% Action and b0 thresholding. Returns set of likely controls under model.
function likely_ctrls = action_and_b0_thresh_ctrls(x, us, gIdx, goals, delta, b0)
    likely_ctrls = {};
    for i=1:length(us)
        u = us{i};
        % Threshold by: P(u | x, g)*b0(g) >= delta
        pu = pugiveng(x, u, gIdx, us, goals)*b0(gIdx);
        if pu >= delta
            likely_ctrls{end+1} = u;
        end
    end
end

%% Belief-based thresholding. Returns set of likely controls under model.
function likely_ctrls = belief_thresh_ctrls(x, us, goals, delta, belief)
    likely_ctrls = {};
    for i=1:length(us)
        u = us{i};
        % Threshold by: \sum_g P(u | x, g)*belief(g) >= delta
        pu = pugiveng(x, u, 1, us, goals)*belief(1) + ...
                pugiveng(x, u, 2, us, goals)*belief(2);
        if pu >= delta
            likely_ctrls{end+1} = u;
        end
    end
end

%% Computes belief update:
%       b1(g) = P(g | x0, u0) = P(u0 | x0; g)b0(g)/\sum_g P(u0 | x0; g)b0(g)
function b1 = next_belief(b0, x0, u0, us, goals)
    normalizer = pugiveng(x0, u0, 1, us, goals)*b0(1) + ...
                    pugiveng(x0, u0, 2, us, goals)*b0(2);
    
    b1 = [0; 0];
        
    b1(1) =  pugiveng(x0, u0, 1, us, goals)*b0(1);
    b1(1) = b1(1)/normalizer;
    
    b1(2) =  pugiveng(x0, u0, 2, us, goals)*b0(2);
    b1(2) = b1(2)/normalizer;
end

%% Computes action model:
%       P(u | x, g) = e^{-||xnext - g||}/\sum_u e^{-||xnext - g||}
function pu = pugiveng(x, u, gIdx, us, goals)
    [xnext, ~] = dynamics(x, u);

    dist_to_goal = norm(goals{gIdx} - xnext);
    
    pu_numerator = exp(-dist_to_goal);
    
    normalizer = 0.;
    for i=1:length(us)
        u = us{i}; 
        [xnext, ~] = dynamics(x, u);
        dist_to_goal = norm(goals{gIdx} - xnext);
        normalizer = normalizer + exp(-dist_to_goal);
    end
    
    pu = pu_numerator / normalizer;
end

%% Computes dynamics.
function [xnext, isvalid] = dynamics(x,u)
    % note: state x = [yidx,xidx] into the 2D prediction matrix
    
    if u == 1 % right
        xnext = [x(1); x(2)+1];
    elseif u == 2 % upright
        xnext = [x(1)-1; x(2)+1];
    elseif u == 3 % up
        xnext = [x(1)-1; x(2)];
    elseif u == 4 % upleft
        xnext = [x(1)-1; x(2)-1];
    elseif u == 5 % left
        xnext = [x(1); x(2)-1];
    elseif u == 6 % leftdown
        xnext = [x(1)+1; x(2)-1];
    elseif u == 7 % down
        xnext = [x(1)+1; x(2)];
    elseif u == 8 % rightdown
        xnext = [x(1)+1; x(2)+1];
    else
        error("Invalid control in dynamics.");
    end
    isvalid = true;
end

%% Plots controls.
function fh = plot_ctrls(x, us, b0, gIdx, goals, fignum)

    fh = figure(fignum);
    hold on;
    for j=1:length(us)
        u = us{j};
        pu = pugiveng(x, u, gIdx, us, goals);
        [xnext, ~] = dynamics(x, u);
        %c = min(1,pu*10);
        c = 0.1; %abs(0.92-c);
        q = quiver(x(1), x(2), x(1)-xnext(1), x(2)-xnext(2), 'Color', [c,c,c]);
        t = text(xnext(1), xnext(2), num2str(pu, 3));
    end
    xlim([1,6]);
    ylim([1,6]);
    scatter(x(1), x(2), 20, 'k');
end

% %% Computes dynamics.
% function [xnext, isvalid] = dynamics(x,u)
%     if u == 1 % right
%         xnext = [x(1)+1; x(2)];
%     elseif u == 2 % upright
%         xnext = [x(1)+1; x(2)+1];
%     elseif u == 3 % up
%         xnext = [x(1); x(2)+1];
%     elseif u == 4 % upleft
%         xnext = [x(1)-1; x(2)+1];
%     elseif u == 5 % left
%         xnext = [x(1)-1; x(2)];
%     elseif u == 6 % leftdown
%         xnext = [x(1)-1; x(2)-1];
%     elseif u == 7 % down
%         xnext = [x(1); x(2)-1];
%     elseif u == 8 % rightdown
%         xnext = [x(1)+1; x(2)-1];
%     else
%         error("Invalid control in dynamics.");
%     end
%     isvalid = true;
% end