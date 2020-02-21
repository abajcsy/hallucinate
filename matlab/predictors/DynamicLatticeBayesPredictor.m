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
        
        states      % (cell arr) indicies of all states in grid
        controls    % (cell arr) all controls
        
        truncG1
        truncG2
        
        rows
        cols
    end
    
    methods
        function obj = DynamicLatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
                gridMin, gridMax, r, hmmParam)
            obj.prior = containers.Map([1:length(prior)], prior);
            obj.goals = goals;
            obj.sigma1 = sigma1;
            obj.sigma2 = sigma2;
            
            obj.gridMin = gridMin;
            obj.gridMax = gridMax;
            
            obj.r = r;
            obj.hmmParam = hmmParam;
                        
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
                        
            pd1 = makedist('Normal','mu',0,'sigma',obj.sigma1);
            pd2 = makedist('Normal','mu',0,'sigma',obj.sigma2);
            obj.truncG1 = truncate(pd1, -pi, pi);
            obj.truncG2 = truncate(pd2, -pi, pi);
            
            % Enumerate all the controls
            %   UP_RIGHT = 1
            %   RIGHT = 2
            %   DOWN_RIGHT = 3
            %   DOWN_LEFT = 4
            %   LEFT = 5
            %   UP_LEFT = 6
            obj.controls = [1,2,3,4,5,6]; 
        end
        
        %% Mega prediction loop for up to timestep H.
        function preds = predict(obj, x0, H)
            
            % Convert into (i,j) index
            [i0, j0] = obj.realToSim(x0);
            
            % Store all the predictions.
            preds = cell([1,H+1]);
            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            preds(:) = {zeros(obj.rows, obj.cols)};
            
            % At current timestep, the measured state has 
            % probability = 1, zeros elsewhere.
            preds{1}(i0, j0) = 1;

            % Make a list that stores the current sequence of all goal
            % values and their corresponding probabilities.
            goalSeq = {[1, obj.prior(1)], [2, obj.prior(2)]};
        
            for t=2:H+1
                % Update the sequence of betas and probabilities for next timestep. 
                [goalSeq, pGoals] = obj.updateGoalSeq(goalSeq);
                fprintf('Predicting for t=%d of T=%d\n', t, H);
                for scurr = obj.states
                    for u = obj.controls
                        % Unpack the [i,j] coords.
                        s = scurr{1};

                        % Optimization!
                        if preds{t-1}(s(1), s(2)) == 0
                            continue;
                        end

                        % now we are cookin' with gas!
                        [snext, isValid] = obj.dynamics(s, u);

                        % if u can take us to a valid state in the
                        % world
                        if isValid
                            for g=1:length(obj.goals)
                                pug = obj.Pu_given_x_g(u, s, g);
                                pgt = pGoals(g);

                                % P(u|x,gt) * P(gt) * P(x)
                                preds{t}(snext(1), snext(2)) = ...
                                    preds{t}(snext(1), snext(2)) + ...
                                    pgt * pug *  preds{t-1}(s(1), s(2));
                            end
                        end
                    end
                end
            end
        end
        
        %% Compute the next betas in the sequence and corresponding probabilities.
        function [newGoalSeq, pGoals] = updateGoalSeq(obj, goalSeq)
            % Store the current probability of goal given all the entire
            % history: 
            %   P(gt) = \sum_{gt-1} ... \sum_{g1}\sum_{g0} P(g0)P(g1 | g0)...P(gt | gt-1)
            % for each gt \in {1,2,..numGoals}
            numGoals = length(obj.goals);
            pGoals = zeros(1, numGoals);
            
            newGoalSeq = {}; 
            for tuple = goalSeq
                % Get the current goal and its probability from the
                % sequence.
                currGoal = tuple{1}(1);
                pCurr = tuple{1}(2);
                
                % Add all possibilities of switching the goals. Compute
                % the probability of switching by querying HMM model. 
                for nextGoal=1:numGoals
                    pNextGivenCurr = obj.hmm(currGoal, nextGoal);
                    pNext = pNextGivenCurr*pCurr;
                    % Update sequence of goals.
                    newGoalSeq{end+1} = [nextGoal, pNext];
                    % Update total probability of each goal.
                    pGoals(nextGoal) = pGoals(nextGoal) + pNext;
                end
            end
        end
        
        %% Encodes HMM dynamic goal dynamics. 
        %  Inputs are the values of the hidden parameter at the next
        %  timestep and the current timestep. The HMM model says that the
        %  value of the parameter at the next time will be the same as the
        %  currrent timstep with probability "hmmParam" and will be
        %  different with probability "1 - hmmParam":
        %       next = {curr with probability (hmmParam + (1-hmmParam)/N)
        %              {next with probability (1-hmmParam)/N
        function prob = hmm(obj, curr, next)
            N = length(obj.goals);
            if curr == next
                prob = obj.hmmParam + (1 - obj.hmmParam)/N;
            else
                prob = (1 - obj.hmmParam)/N;
            end
        end
        
        %% Helper function. 
        function plot_Pu_given_x(obj)
            
            gString = createGrid(obj.gridMin, obj.gridMax, obj.gridDims);
            
            for u = obj.controls
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
        function prob = Pu_given_x_g(obj, u, s0, goal)
            % Get the lower and upper bounds to integrate the Gaussian 
            % PDF over.
            bounds = obj.uToThetaBounds(u);
            [x, y] = obj.simToReal(s0);        
            
            if goal == 1
                % Compute optimal control (i.e. mean of Gaussian) 
                g1 = obj.goals{1};
                mu1 = atan2(g1(2) - y, g1(1) - x); 
                
                % Find the integration bounds. 
                bound1 = wrapToPi(mu1 - bounds(1));
                bound2 = wrapToPi(mu1 - bounds(2));
                
                % Integrate on bounds. 
                p = cdf(obj.truncG1, [bound1, bound2]);
                prob = abs(p(2) - p(1));
            elseif goal == 2
                % Compute optimal control (i.e. mean of Gaussian) 
                g2 = obj.goals{2};
                mu2 = atan2(g2(2) - y, g2(1) - x); 
                
                % Find the integration bounds. 
                bound1 = wrapToPi(mu2 - bounds(1));
                bound2 = wrapToPi(mu2 - bounds(2));
                
                % Integrate on bounds. 
                p = cdf(obj.truncG2, [bound1, bound2]);
                prob = abs(p(2) - p(1));
            else
                error("In PuGivenGoal(): goal is invalid: %d\n", goal);
            end
            
            
            if abs(bound1 - bound2) > pi
                prob = 1 - prob;
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
        
        %% Converts from fake discrete controls into lower and upper theta
        function bounds = uToThetaBounds(obj, u)
            if u == 1 % UP_RIGHT
                bounds = [pi/6, pi/2];
            elseif u == 2 % RIGHT
                bounds = [-pi/6, pi/6];
            elseif u == 3 % DOWN_RIGHT
                bounds = [-pi/2, -pi/6];
            elseif u == 4 % DOWN_LEFT
                bounds = [-(5*pi)/6, -pi/2];
            elseif u == 5 % LEFT
                bounds = [-(5*pi)/6, (5*pi)/6];
            elseif u == 6 % UP_LEFT
                bounds = [pi/2, (5*pi)/6];
            end
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
