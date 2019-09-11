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
        dxEven
        dxOdd
        dy
        
        states      % (cell arr) indicies of all states in grid
        controls    % (cell arr) all controls
        
        truncG1
        truncG2
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
            obj.dxEven = obj.r;
            obj.dxOdd = obj.r/2;
            obj.dy = (obj.r*sqrt(3))/2.;
                        
            % Enumerate all the state indicies
            obj.states = {};
             for y = obj.gridMin(2):obj.dy:obj.gridMax(2)
                i = round((y - obj.gridMin(2))/obj.dy);
                if mod(i,2) == 0
                    for x = obj.gridMin(1):obj.dxEven:obj.gridMax(1)
                        obj.states{end+1} = obj.realToSim([x,y]);
                    end
                else
                    for x = obj.gridMin(1):obj.dxOdd:obj.gridMax(1)
                        obj.states{end+1} = obj.realToSim([x,y]);
                    end
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
        function fullPreds = predict(obj, x0, H)
            
            % Convert into (i,j) index
            coor = obj.realToSim(x0);
            
            % Make a map with keys being goals and values being
            % predictions.
            twoGoalPreds = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            for g=1:length(obj.goals)
                predsG = cell([1,H+1]);
                predsG(:) = {zeros(obj.gridDims)};
                
            	% At current timestep, the measured state has 
                % probability = 1, zeros elsewhere.
                predsG{1}(coor(1), coor(2)) = 1;
                
                twoGoalPreds(num2str(g)) = predsG;
            end
            
            % First compute all the predictions for goal 1.
            for g = 1:length(obj.goals)
                preds = twoGoalPreds(num2str(g));
                for t=2:H+1
                    fprintf('Predicting for g=%d for t=%d\n', g, t);
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
                                pug = obj.Pu_given_x_g(u, s, g);

                                % P(u|x,g) * P(g) * P(x)
                                preds{t}(snext(1), snext(2)) = ...
                                    preds{t}(snext(1), snext(2)) + ...
                                    pug *  preds{t-1}(s(1), s(2));
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
            fullPreds(:) = {zeros(obj.gridDims)};
            fullPreds{1}(coor(1), coor(2)) = 1;
            for t=2:H+1
                fullPreds{t} = predsG1{t}*obj.prior(1) + predsG2{t}*obj.prior(2);
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
            x = obj.simToReal(s0);        
            
            if goal == 1
                % Compute optimal control (i.e. mean of Gaussian) 
                g1 = obj.goals{1};
                mu1 = atan2(g1(2) - x(2), g1(1) - x(1)); 
                
                % Find the integration bounds. 
                bound1 = wrapToPi(mu1 - bounds(1));
                bound2 = wrapToPi(mu1 - bounds(2));
                
                % Integrate on bounds. 
                p = cdf(obj.truncG1, [bound1, bound2]);
                prob = abs(p(2) - p(1));
            elseif goal == 2
                % Compute optimal control (i.e. mean of Gaussian) 
                g2 = obj.goals{2};
                mu2 = atan2(g2(2) - x(2), g2(1) - x(1)); 
                
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
        function  coor = realToSim(obj, x)
            i = round((x(2) - obj.gridMin(2))/(obj.r*sqrt(3)/2.));
            j = round((x(1) - obj.gridMin(1))/obj.r);
            
            if mod(i,2) ~= 0 % if odd, need to shift by r/2
                j = j + obj.r/2.;
            end
            
            coor = [i,j];
        end
        
        %% Converts from coordinate (i,j) to real state (m)
        function real = simToReal(obj, coor)
            y = coor(1)*(obj.r * sqrt(3)/2.) + obj.gridMin(2);
            
            if mod(coor(1), 2) == 0
                x = coor(2)*obj.r + obj.gridMin(1);
            else
                x = (coor(2) - obj.r/2.)*obj.r + obj.gridMin(1);
            end
            
            real = [x,y];
        end
        
        %% Converts from fake discrete controls into lower and upper theta
        function bounds = uToThetaBounds(obj, u)
            if u == 1 % UP_RIGHT
                bounds = [pi/6, pi/2];
            elseif u == 2 % RIGHT
                bounds = [-pi/6, pi/2];
            elseif u == 3 % DOWN_RIGHT
                bounds = [-pi/6, -pi/2];
            elseif u == 4 % DOWN_LEFT
                bounds = [-pi/2, -(5*pi)/6];
            elseif u == 5 % LEFT
                bounds = [(5*pi)/6, -(5*pi)/6];
            elseif u == 6 % UP_LEFT
                bounds = [pi/2, (5*pi)/6];
            end
        end
        
        %% Dynamics function gives next state we can get to given current state.
        function [snext, isValid] = dynamics(obj, s0, u)
            isValid = true;
            
            % Note: s = [i, j] 
            % where i represents "rows" (y), j represents "columns" (x)
            
            if u == 1 % UP_RIGHT
                snext = [s0(1)+1; s0(2)+1];
            elseif u == 2 % RIGHT
                snext = [s0(1); s0(2)+1];
            elseif u == 3 % DOWN_RIGHT
                snext = [s0(1)-1; s0(2)];
            elseif u == 4 % DOWN_LEFT
                snext = [s0(1)-1; s0(2)-1];
            elseif u == 5 % LEFT
                snext = [s0(1); s0(2)-1];
            elseif u == 6 % UP_LEFT
                snext = [s0(1); s0(2)-1];
            end
            
            if snext(1) < 1 || snext(1) > obj.gridDims(1) || ...
                    snext(2) < 1 || snext(2) > obj.gridDims(2)  
                snext = s0;
                isValid = false;
            end
        end
        
    end
end

