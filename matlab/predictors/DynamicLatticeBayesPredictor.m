classdef DynamicLatticeBayesPredictor
    %DYNAMICLATTICEBAYESPREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma1      % sigma in Gaussian observation model
        sigma2      
        
        goals       % (array) possible goals
        prior       % (map) prior over goals -- goal number is key, probability is value
        
        gridMin
        gridMax
        
        r           % Grid cell size.
        hmmParam    % Probability of unknown parameter's value staying the same.
        Delta       % Discrete distribution. Probability of drawing goal \Delta{1,2}. 
        
        states      % (cell arr) indicies of all states in grid
        us          % (cell arr) real-world angles of each control
        usIdxs      % (cell arr) index of each control
        ubounds     % (cell arr) integration bounds for each control.
        
        %truncG1
        %truncG2
        truncpd     % truncated zero-mean gaussian.
        
        rows
        cols
    end
    
    methods
        function obj = DynamicLatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
                gridMin, gridMax, r, hmmParam, Delta)
            obj.prior = containers.Map([1:length(prior)], prior);
            obj.goals = goals;
            obj.sigma1 = sigma1;
            obj.sigma2 = sigma2;
            
            obj.gridMin = gridMin;
            obj.gridMax = gridMax;
            
            obj.r = r;
            obj.hmmParam = hmmParam;
            obj.Delta = Delta;
                        
            % Enumerate all the state indicies
            obj.states = {};
            
            % Number of cells along y (i.e. number of rows).
            obj.rows = round((obj.gridMax(2) - obj.gridMin(2)) / ...
                (0.5 * obj.r * sqrt(3)));
            
            % Number of cells along x (i.e. number of columns).
            obj.cols = round((obj.gridMax(1) - obj.gridMin(1)) / ...
                obj.r); 
            
            for i = 1:obj.rows
                for j = 1:obj.cols
                    obj.states{end + 1} = [i, j];
                end
            end
            
            pd = makedist('Normal','mu',0,'sigma',obj.sigma1);
            obj.truncpd = truncate(pd, -pi, pi);
            
            % Enumerate all the controls for a lattice predictor.
            %   UP_RIGHT = 1
            %   RIGHT = 2
            %   DOWN_RIGHT = 3
            %   DOWN_LEFT = 4
            %   LEFT = 5
            %   UP_LEFT = 6
            obj.usIdxs = [1, 2, 3, 4, 5, 6]; 
            obj.us = [pi/3, 0, -pi/3, -(2*pi)/3, -pi, (2*pi)/3];
            obj.ubounds = obj.genUBounds();
            
        end
        
        %% Mega prediction loop for up to timestep H.
        function fullPreds = predict(obj, x0, H)
            
            % Convert into (i,j) index
            [i0, j0] = obj.realToSim(x0);
            
            % At current timestep, the measured state has 
            % probability = 1, zeros elsewhere.
            px0 = zeros(obj.rows, obj.cols);
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
                    
            for t=2:H+1                
                % Update the goal sequence for current t.
                goalSeq = obj.updateGoalSeq(goalSeq, t);
                
                % Get the set of all goal sequences 
                %   (e.g. {{1,1}, {2,1}, ... {2,2}}
                gtseq = keys(goalSeq(num2str(t)));
                
                fprintf('Predicting for t=%d of T=%d\n', t, H);
                
                % Get the map with all the probabilities at the prior time.
                p_xtm1 = predSeq(num2str(t-1));

                % Create map for current time. 
                p_xt = containers.Map('KeyType', 'char', 'ValueType', 'any');
                
                for gs_cell = gtseq
                    curr_gs = str2num(gs_cell{1});
                    prev_gs = curr_gs(1:end-1);

                    % catch corner case when we are looking back to x0.
                    if isempty(prev_gs)
                        prev_gs = [1];
                    end

                    % NOTE: THIS IS HARDCODED FOR NOW FOR 2 GOALS
                    p_xt_prev = p_xtm1(num2str(prev_gs));

                    % Get which goal was added this time .
                    gIdx = curr_gs(end);
                    
                    % Do the prediction!
                    preds_curr = obj.updateStateDist(p_xt_prev, gIdx);
                        
                    % Store the prediction!
                    p_xt(num2str(curr_gs)) = preds_curr;
                end
                % Add the predictions for all mode sequences to current time.
                predSeq(num2str(t)) = p_xt;
            end
            
            % Combine all predictions over all sequences
            fullPreds = cell([1,H+1]);
            fullPreds(:) = {zeros(obj.rows, obj.cols)};
            fullPreds{1}(i0, j0) = 1;
            for t=2:H+1
                % Get goal sequence at this time.
               curr_seqs = goalSeq(num2str(t));
               curr_preds = predSeq(num2str(t));
               gtseq = keys(curr_seqs);
               for gs_cell = gtseq
                   g_seq = str2num(gs_cell{1});
                   p_g_seqt = curr_seqs(num2str(g_seq));
                   p_x_g_seqt = curr_preds(num2str(g_seq));
                   fullPreds{t} = fullPreds{t} + p_x_g_seqt * p_g_seqt;
               end
            end
        end
        
        
        %% Updates state distribution given a specific goal and 
        % the prior predictions. 
        function preds_curr = updateStateDist(obj, preds_prev, curr_gIdx)
            preds_curr = zeros(obj.rows, obj.cols);
            for scurr = obj.states
                for uid = obj.usIdxs
                    % Unpack the [i,j] coords.
                    s = scurr{1};

                    % Optimization!
                    if preds_prev(s(1), s(2)) == 0
                        continue;
                    end

                    % now we are cookin' with gas!
                    [snext, isValid] = obj.dynamics(s, uid);

                    % if u can take us to a valid state in the
                    % world
                    if isValid
                        ureal = obj.us(uid);
                        pug = obj.Pu_given_x_g(ureal, s, curr_gIdx);

                        % P(x_t+1 | x_t, g)
                        preds_curr(snext(1), snext(2)) = ...
                            preds_curr(snext(1), snext(2)) + ...
                            pug *  preds_prev(s(1), s(2));
                    end
                end
            end
        end
        
        %% Updates the goal sequence. 
        function goalSeq = updateGoalSeq(obj, goalSeq, t)
             if isempty(goalSeq)
                % Setup t = 2 priors. 
                g1seq = containers.Map('KeyType', 'char', 'ValueType', 'any');
                g1seq(num2str([1])) = obj.prior(1);
                g1seq(num2str([2])) = obj.prior(2);
                goalSeq(num2str(t)) = g1seq;
                return; 
            end

            g_prev_map = goalSeq(num2str(t-1));
            g_prev_keys = keys(g_prev_map);
            g_curr_map = ...
                containers.Map('KeyType', 'char', 'ValueType', 'any');
            for gs_cell = g_prev_keys
                prev_gs = str2num(gs_cell{1});
                for currGoal=1:length(obj.goals)
                    % Get the previous goal we had. 
                    prevGoal = prev_gs(end);
                    % Get the probability of getting to that previous sequence. 
                    pPrev = g_prev_map(num2str(prev_gs));

                    % Setup the new goal sequence.
                    curr_gs = prev_gs;
                    curr_gs(end+1) = currGoal;

                    % Get P(curr | prev)
                    pCurrGivenPrev = obj.hmm(currGoal, prevGoal);

                    % Get P(curr) = P(curr | prev) * P(prev)
                    pCurr = pCurrGivenPrev * pPrev;

                    % Set the probability of this new expanded goal sequence. 
                    g_curr_map(num2str(curr_gs)) = pCurr;
                end
            end
            % Add all the new sequences to the corresponding time in map. 
            goalSeq(num2str(t)) = g_curr_map;
        end

        %% Encodes HMM dynamic goal dynamics. 
        %  Inputs are the values of the hidden parameter at the curr
        %  timestep and the previous timestep. The HMM model says that the
        %  value of the parameter at the curr time will be the same as the
        %  previous timstep with probability "hmmParam" and will be
        %  different with probability "1 - hmmParam":
        % 
        %       curr = {prev with probability Delta*(1-hmmParam)
        %              {curr with probability hmmParam + Delta*(1-hmmParam)
        % 
        function prob = hmm(obj, curr, prev)
            if curr == prev
                prob = obj.hmmParam + obj.Delta(curr)*(1-obj.hmmParam);
            else
                prob = obj.Delta(curr)*(1-obj.hmmParam);
            end
        end
        
        %% Helper function. 
        function plot_Pu_given_x(obj)
            
            gString = createGrid(obj.gridMin, obj.gridMax, obj.gridDims);
            
            for u = obj.us
                grid = zeros(obj.gridDims);
                for scurr = obj.states
                    s = scurr{1};
                    pug1 = obj.Pu_given_x_g(u, s, 1);
                    pug2 = obj.Pu_given_x_g(u, s, 2);
                    grid(s(1), s(2)) = pug1*obj.prior(1) + pug2*obj.prior(2);
                end
                pcolor(gString.xs{2}, gString.xs{1}, grid);
                title(strcat('P(u=', num2str(u), '|x)'));
                colorbar
                caxis([0,1])
            end

        end
        
        %% Compute the probability of a specific action given a state and beta value.
        %  Gaussian observation model:
        %
        %       P(u | x0; g1) = N(atan2(g1(y) - y, g1(x) - x), sigma_1^2)
        %
        %  and 
        % 
        %       P(u | x0; g2) = N(atan2(g2(y) - y, g2(x) - x), sigma_2^2)
        %
        function prob = Pu_given_x_g(obj, u, s0, goalIdx)
            % Get the lower and upper bounds to integrate the Gaussian 
            % PDF over.
            [x, y] = obj.simToReal(s0);    
            uopt = atan2(obj.goals{goalIdx}(2)- y, obj.goals{goalIdx}(1) - x);
           
            % minimum angular distance between current control (u) and uopt
            diff = abs(angdiff(u, uopt));

            % corner case.
            if length(obj.us) == 1
                prob = 1;
                return
            end

            % note because of numerics: 
            % sometimes we get controls that are just close to zero but are negative
            zero_tol = -1e-7; 
            % find indicies of all controls in the positive [0, pi] range.
            pos_idxs = find(obj.us >= zero_tol);
            pos_idxs(end+1) = find(obj.us == -pi); %include -pi since its = pi

            % find the control bounds.
            for i=pos_idxs
                bounds = obj.ubounds{i};
                low_bound = bounds(1);
                up_bound = bounds(2);

                if low_bound < 0 && diff <= up_bound
                    % catch corner case around 0.
                    p = cdf(obj.truncpd, [0, up_bound]);
                    prob = abs(p(2) - p(1)) * 2;
                elseif up_bound < 0 && diff >= low_bound
                    % considering an action that falls into the bounds around -pi
                    p = cdf(obj.truncpd, [low_bound, pi]);
                    prob = abs(p(2) - p(1)) * 2;
                elseif diff <= up_bound && diff >= low_bound
                    % normal integration over bounds.
                    p = cdf(obj.truncpd, [low_bound, up_bound]);
                    prob = abs(p(2) - p(1));
                end
            end
        end
        
        %% Converts from real state (m) to coordinate (i,j)
        function  [i, j] = realToSim(obj, x)
            i = round((2 * (x(2) - obj.gridMin(2))) / (obj.r * sqrt(3))) + 1;
            
            if mod(i, 2) == 1
                % Case where i is odd.
                j = round((x(1) - obj.gridMin(1)) / obj.r) + 1;
            else
                % Case where i is even.
                j = round((x(1) - obj.gridMin(1) - 0.5 * obj.r) / obj.r) + 2;
            end
        end
        
        %% Converts from coordinate (i,j) to real state (m)
        function [x, y] = simToReal(obj, coor)
            i = coor(1);
            j = coor(2);
            
            y = obj.gridMin(2) + (i - 1) * (obj.r * 0.5 * sqrt(3));
            
            if mod(i, 2) == 1
                x = obj.gridMin(1) + (j - 1) * obj.r;
            else
                x = obj.gridMin(1) + 0.5 * obj.r + obj.r * (j - 2);
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
        
        %% Dynamics function gives next state we can get to given current state.
        function [snext, isValid] = dynamics(obj, s0, u)
            isValid = true;
            
            % Note: s = [i, j] 
            % where i represents "rows" (y), j represents "columns" (x)
            i = s0(1);
            j = s0(2);
            
            if u == 1 % UP_RIGHT
                if mod(i, 2) == 0
                    snext = [i + 1, j];
                else
                    snext = [i + 1, j + 1];
                end
            elseif u == 2 % RIGHT
                snext = [i, j + 1];
            elseif u == 3 % DOWN_RIGHT
                if mod(i, 2) == 0
                    snext = [i - 1, j];
                else
                    snext = [i - 1, j + 1];
                end
            elseif u == 4 % DOWN_LEFT
                if mod(i, 2) == 0
                    snext = [i - 1, j - 1];
                else
                    snext = [i - 1, j];
                end
            elseif u == 5 % LEFT
                snext = [i, j - 1];
            elseif u == 6 % UP_LEFT
                if mod(i, 2) == 0
                    snext = [i + 1, j - 1];
                else
                    snext = [i + 1, j];
                end
            end
            
            if snext(1) < 1 || snext(1) > obj.rows || ...
                    snext(2) < 1 || snext(2) > obj.cols  
                snext = s0;
                isValid = false;
            end
        end
        
        %% Gets the meshgrid that represents the real-world coords.
        function [X, Y] = getLatticeMeshgrid(obj)
            X = zeros(obj.rows, obj.cols);
            Y = zeros(obj.rows, obj.cols);
            for i = 1:obj.rows
                for j = 1:obj.cols
                    [x, y] = obj.simToReal([i, j]);
                    X(i, j) = x;
                    Y(i, j) = y;
                end
            end
        end
        
    end
end
