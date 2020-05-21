classdef ThresholdedLatticeBayesPredictor
    %THRESHOLDEDLATTICEBAYESPREDICTOR Summary of this class goes here
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
        
        %truncG1
        %truncG2
        truncpd     % truncated zero-mean gaussian.
        
        rows
        cols
        
        uThresh
    end
    
    methods
        function obj = ThresholdedLatticeBayesPredictor(prior, ...
                goals, sigma1, sigma2, ...
                gridMin, gridMax, r, uThresh)
            obj.prior = containers.Map([1:length(prior)], prior);
            obj.goals = goals;
            obj.sigma1 = sigma1;
            obj.sigma2 = sigma2;
            
            obj.gridMin = gridMin;
            obj.gridMax = gridMax;
            
            obj.r = r;
            obj.uThresh = uThresh;
                        
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
                        
%             pd1 = makedist('Normal','mu',0,'sigma',obj.sigma1);
%             pd2 = makedist('Normal','mu',0,'sigma',obj.sigma2);
%             obj.truncG1 = truncate(pd1, -pi, pi);
%             obj.truncG2 = truncate(pd2, -pi, pi);

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
            
            % Initialize empty prediction grids forward in time.
            fullPreds = cell([1,H+1]);
            fullPreds(:) = {zeros(obj.rows, obj.cols)};
            % At current timestep, the measured state has 
            % probability = 1, zeros elsewhere.
            fullPreds{1}(i0, j0) = 1;
            
            for t=2:H+1
                for scurr = obj.states
                    for uid = obj.usIdxs
                        % Unpack the [i,j] coords.
                        s = scurr{1};

                        % Optimization!
                        if fullPreds{t-1}(s(1), s(2)) == 0
                            continue;
                        end

                        % now we are cookin' with gas!
                        [snext, isValid] = obj.dynamics(s, uid);

                        % if u can take us to a valid state in the
                        % world
                        if isValid
                            ureal = obj.us(uid);
                            
                            % Compute:
                            %   P(u | x) = \sum_g P(u | x, g)P(g)
                            pu = 0.0;
                            for gidx=1:length(obj.goals)
                               pu = pu + ...
                                   obj.Pu_given_x_g(ureal, s, gidx)*obj.prior(gidx);
                            end
                            
%                             % Apply thresholding to controls!
%                             if pu >= obj.uThresh
%                                 pu = 1.0;
%                             else
%                                 pu = 0.0;
%                             end

                            % Compute:
                            %   P(x' | x) = \sum_u P(x' | x, u) * P(u|x) * P(x)
                            fullPreds{t}(snext(1), snext(2)) = ...
                                fullPreds{t}(snext(1), snext(2)) + ...
                                pu *  fullPreds{t-1}(s(1), s(2));
                        end
                    end
                end
%                 % Make this a binary occupancy map.
%                 fullPreds{t}(fullPreds{t} > 0) = 1;
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

