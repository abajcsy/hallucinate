% Enumerate all the controls for a lattice predictor.
%   UP_RIGHT = 1
%   RIGHT = 2
%   DOWN_RIGHT = 3
%   DOWN_LEFT = 4
%   LEFT = 5
%   UP_LEFT = 6
usIdxs = [1, 2, 3, 4, 5, 6]; 
us = [pi/3, 0, -pi/3, -(2*pi)/3, -pi, (2*pi)/3];
ubounds = genUBounds();

% At current timestep, the measured state has 
% probability = 1, zeros elsewhere.
i0 = 47;
j0 = 41;
rows = 92;
cols = 80;
px0 = zeros(rows, cols);
px0(i0, j0) = 1;

% Data structure for *predictions*.
% Map with Keys = times
%          Value = Maps with Keys = goal sequence
%                            Values = predictions
% t = 1 ---> {1} --> p(x0)
%       ---> {2} --> p(x0)
% t = 2 ---> {1} --> p(x1|g1=1)
%            {2} --> p(x1|g1=2)
% t = 3 ---> {1,1} --> p(x2|g2=1,g1=1)
%            {2,1} --> p(x2|g2=1,g1=2)
%            {1,2} --> p(x2|g2=2,g1=1)
%            {2,2} --> p(x2|g2=2,g1=2)
% t = 4 ---> ....
predSeq = containers.Map('KeyType', 'char', 'ValueType', 'any');
% Setup t = 1 state distributions. 
px0dists = containers.Map('KeyType', 'char', 'ValueType', 'any');
px0dists(num2str([1])) = px0;
px0dists(num2str([2])) = px0;
predSeq(num2str(1)) = px0dists;

% Data structure for *goal sequence*. 
% Map with Keys = times
%          Value = Maps with Keys = goal sequence
%                            Values = probability of goal seq
% t = 2 ---> {1} --> p(g1=1)
%            {2} --> p(g1=2)
% t = 3 ---> {1,1} --> p(g2=1 | g1=1)*P(g1=1)
%            {2,1} --> p(g2=1 | g1=2)*P(g1=2)
%            {1,2} --> p(g2=2 | g1=1)*P(g1=1)
%            {2,2} --> p(g2=2 | g1=2)*P(g1=2)
% t = 4 ---> ....
goalSeq = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Setup P(g0 = 1) and P(g0 = 2)
prior = [0.5, 0.5];
H = 3;
for t=2:H
    % Update the goal sequence for current t.
    goalSeq = updateGoalSeq(goalSeq, t, prior);
    % Get the set of all goal sequences 
    %   (e.g. {{1,1}, {2,1}, ... {2,2}}
    gtseq = keys(goalSeq(num2str(t)));
    
    % Get the map with all the probabilities at the prior time.
    p_xtm1 = predSeq(num2str(t-1));
        
    % Create map for current time. 
    p_xt = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for gs_cell = gtseq
        curr_gs = str2num(gs_cell{1})
        prev_gs = curr_gs(1:end-1);
        
        % catch corner case when we are looking back to x0.
        if isempty(prev_gs)
            prev_gs = [1];
        end
        
        p_xt_prev = p_xtm1(num2str(prev_gs));
        for g=1:2
            curr_pred = zeros(rows, cols);
            
            % --- do the prediction here. --- %
            
            p_xt(num2str(curr_gs)) = curr_pred;
        end
    end
    % Add the predictions for all mode sequences to current time.
    predSeq(num2str(t)) = p_xt;
end

% Combine all predictions over all sequences
fullPreds = cell([1,H+1]);
fullPreds(:) = {zeros(rows, cols)};
fullPreds{1}(i0, j0) = 1;
for t=2:H
    % Get goal sequence at this time.
   curr_seqs = goalSeq(num2str(t));
   curr_preds = predSeq(num2str(t));
   gtseq = keys(curr_seqs);
   for gs_cell = gtseq
       seq = str2num(gs_cell{1});
       pseqt = curr_seqs(num2str(seq));
       pxseqt = curr_preds(num2str(seq));
       fullPreds{t} = fullPreds{t} + pxseqt*pseqt;
   end
end

%% Updates goal seq at current time. 
function goalSeq = updateGoalSeq(goalSeq, t, prior)
    if isempty(goalSeq)
        % Setup t = 2 priors. 
        g1seq = containers.Map('KeyType', 'char', 'ValueType', 'any');
        g1seq(num2str([1])) = prior(1);
        g1seq(num2str([2])) = prior(2);
        goalSeq(num2str(t)) = g1seq;
        return; 
    end
    
    numGoals = 2;
    g_prev_map = goalSeq(num2str(t-1));
    g_prev_keys = keys(g_prev_map);
    g_curr_map = ...
        containers.Map('KeyType', 'char', 'ValueType', 'any');
    for gs_cell = g_prev_keys
        prev_gs = str2num(gs_cell{1});
        for currGoal=1:numGoals
            % Get the previous goal we had. 
            prevGoal = prev_gs(end);
            % Get the probability of getting to that previous sequence. 
            pPrev = g_prev_map(num2str(prev_gs));
            
            % Setup the new goal sequence.
            curr_gs = prev_gs;
            curr_gs(end+1) = currGoal;
            
            % Get P(curr | prev)
            pCurrGivenPrev = hmm(currGoal, prevGoal);
            
            % Get P(curr) = P(curr | prev) * P(prev)
            pCurr = pCurrGivenPrev * pPrev;
            
            % Set the probability of this new expanded goal sequence. 
            g_curr_map(num2str(curr_gs)) = pCurr;
        end
    end
    % Add all the new sequences to the corresponding time in map. 
    goalSeq(num2str(t)) = g_curr_map;
end

%% HMM model (with fixed alpha). 
function p = hmm(curr, prev)
    if curr == prev
        p = 1.0;
    else
        p = 0.0;
    end
end

%% Gets the control bounds for integration. 
function ubounds = genUBounds(obj)
    ubounds{1} = [pi/6, pi/2];              % UP_RIGHT
    ubounds{2} = [-pi/6, pi/6];             % RIGHT
    ubounds{3} = [-pi/2, -pi/6];            % DOWN_RIGHT
    ubounds{4} = [-(5*pi)/6, -pi/2];        % DOWN_LEFT
    ubounds{5} = [(5*pi)/6, -(5*pi)/6];     % LEFT
    ubounds{6} = [pi/2, (5*pi)/6];          % UP_LEFT
end