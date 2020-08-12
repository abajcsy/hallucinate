classdef LatticeBayesPredictor
    %LATTICEBAYESPREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma1      % sigma in Gaussian observation model
        sigma2      
        
        goals       % (array) possible goals
        prior       % (map) prior over goals -- goal number is key, probability is value
        
        gridMin
        gridMax
        
        r
        
        states      % (cell arr) indicies of all states in grid
        us          % (cell arr) real-world angles of each control
        usIdxs      % (cell arr) index of each control
        ubounds     % (cell arr) integration bounds for each control.
        
        truncG1
        truncG2
        %truncpd     % truncated zero-mean gaussian.
        
        rows
        cols
    end
    
    methods
        function obj = LatticeBayesPredictor(prior, goals, sigma1, sigma2, ...
                gridMin, gridMax, r)
            obj.prior = containers.Map([1:length(prior)], prior);
            obj.goals = goals;
            obj.sigma1 = sigma1;
            obj.sigma2 = sigma2;
            
            obj.gridMin = gridMin;
            obj.gridMax = gridMax;
            
            obj.r = r;
                        
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

%             pd = makedist('Normal','mu',0,'sigma',obj.sigma1);
%             obj.truncpd = truncate(pd, -pi, pi);
            
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
            
            % Make a map with keys being goals and values being
            % predictions.
            twoGoalPreds = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            for g=1:length(obj.goals)
                predsG = cell([1,H+1]);
                predsG(:) = {zeros(obj.rows, obj.cols)};
                
            	% At current timestep, the measured state has 
                % probability = 1, zeros elsewhere.
                predsG{1}(i0, j0) = 1;
                
                twoGoalPreds(num2str(g)) = predsG;
            end
            
            % First compute all the predictions for goal 1.
            for g = 1:length(obj.goals)
                preds = twoGoalPreds(num2str(g));
                for t=2:H+1
%                     fprintf('Predicting for g=%d for t=%d\n', g, t);
                    for scurr = obj.states
                        for uid = obj.usIdxs
                            % Unpack the [i,j] coords.
                            s = scurr{1};

                            % Optimization!
                            if preds{t-1}(s(1), s(2)) == 0
                                continue;
                            end

                            % now we are cookin' with gas!
                            [snext, isValid] = obj.dynamics(s, uid);

                            % if u can take us to a valid state in the
                            % world
                            if isValid
                                ureal = obj.us(uid);
                                pug = obj.Pu_given_x_g_normalized(ureal, s, g);

                                % P(x'|x, g) = \sum_u P(x'|x,u) * P(u|x,g) * P(x)
                                preds{t}(snext(1), snext(2)) = ...
                                    preds{t}(snext(1), snext(2)) + ...
                                    pug * preds{t-1}(s(1), s(2));
                            end
                        end
                    end
                end
                twoGoalPreds(num2str(g)) = preds;
            end
            
            % Combine the two predictions for each goal by multiplying with
            % the prior over each goal. 
            predsG1 = twoGoalPreds(num2str(1));
            predsG2 = twoGoalPreds(num2str(2));
            
            fullPreds = cell([1,H+1]);
            fullPreds(:) = {zeros(obj.rows, obj.cols)};
            fullPreds{1}(i0, j0) = 1;
            for t=2:H+1
                fullPreds{t} = predsG1{t}*obj.prior(1) + predsG2{t}*obj.prior(2);
            end
            
            % ---- debugging ----%
%             t=3; % time when to sanity check. 
%             linidxs = find(preds{t} > 0.0); 
%             for i =1:length(linidxs)
%                 lidx = linidxs(i);
%                 [row,col] = ind2sub([obj.rows, obj.cols],lidx);
%                 [x,y] = obj.simToReal([row, col]);
%                 fprintf(strcat("P(x(",num2str(t),")=", ...
%                     num2str(x), ",", num2str(y), ") =", num2str(fullPreds{t}(lidx)),"\n")); 
%             end
            % ---- debugging ----%
        end
        
        %% Helper function. 
        function plot_Pu_given_x(obj)
            
            gString = createGrid(obj.gridMin, obj.gridMax, obj.gridDims);
            
            for u = obj.us
                grid = zeros(obj.gridDims);
                for scurr = obj.states
                    s = scurr{1};
                    %pug1 = obj.Pu_given_x_g(u, s, 1);
                    %pug2 = obj.Pu_given_x_g(u, s, 2);
                    pug1 = obj.Pu_given_x_g_normalized(u, s, 1);
                    pug2 = obj.Pu_given_x_g_normalized(u, s, 2);
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
        
        
        %% Compute the probability of a specific action given a state and beta value.
        %  Gaussian observation model:
        %
        %       P(u | x0; g1) = N(atan2(g1(y) - y, g1(x) - x), sigma_1^2)
        %
        %  and 
        % 
        %       P(u | x0; g2) = N(atan2(g2(y) - y, g2(x) - x), sigma_2^2)
        %
        function prob = Pu_given_x_g_normalized(obj, u, s0, goalIdx)
            % Get the lower and upper bounds to integrate the Gaussian 
            % PDF over.
            [x, y] = obj.simToReal(s0);    
            uopt = atan2(obj.goals{goalIdx}(2)- y, obj.goals{goalIdx}(1) - x);
            
            % minimum angular distance between current control (u) and uopt
            diff = abs(angdiff(u, uopt));

            truncpd = [];
            if goalIdx == 1
                truncpd = obj.truncG1;
            elseif goalIdx == 2
                truncpd = obj.truncG2;
            else
                error("Goal idx is not valid in Pugiveng.");
            end
            
            % Get all the probabilities
            unnormalized_probs = pdf(truncpd, diff);
            
            % normalize.
            norm = zeros(size(u));
            for otheru=obj.us
                otheru_arr = otheru * ones(size(u));
                otherdiff = abs(angdiff(otheru_arr, uopt));
                new_prob = pdf(truncpd, otherdiff);
                norm = norm + new_prob;
            end
            
            prob = unnormalized_probs ./ norm;
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

