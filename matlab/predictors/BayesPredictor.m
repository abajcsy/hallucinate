classdef BayesPredictor
    %BAYESPREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sigma1      % sigma in Gaussian observation model
        sigma2      
        
        goals       % (array) possible goals
        prior       % (map) prior over goals -- goal number is key, probability is value
        
        gridMin
        gridMax
        gridDims    % (array) num grid cells in each dim (i.e. height (rows) and width (cols))
        
        dx          % (float) resolution in x and y
        dy
        
        states      % (cell arr) indicies of all states in grid
        controls    % (cell arr) all controls
        
        cdfs        % (Map) indexed by string that has 'goalNum, i, j'
    end
    
    methods
        function obj = BayesPredictor(prior, goals, sigma1, sigma2, gridMin, gridMax, gridDims)
            obj.prior = containers.Map([1:length(prior)], prior);
            obj.cdfs = containers.Map();
            obj.goals = goals;
            obj.sigma1 = sigma1;
            obj.sigma2 = sigma2;
            
            obj.gridDims = gridDims;
            obj.gridMin = gridMin;
            obj.gridMax = gridMax;
            
            obj.dx = (obj.gridMax(1) - obj.gridMin(1))/(obj.gridDims(2)-1);
            obj.dy = (obj.gridMax(2) - obj.gridMin(2))/(obj.gridDims(1)-1);
                        
            % Enumerate all the state indicies
            obj.states = {};
            for i=1:obj.gridDims(1)
                for j=1:obj.gridDims(2)
                    obj.states{end+1} = [i, j];
                    x = obj.simToReal([i, j]);
                    
                    % Compute CDF for each state for goal 1
                    g1 = obj.goals{1};
                    mu1 = atan2(g1(2)- x(2), g1(1) - x(1)); 
                    pd1 = makedist('Normal','mu',mu1,'sigma',obj.sigma1);
                    trunc1 = truncate(pd1, -pi, pi);
                    
                    key1 = strcat(num2str(1), ',', num2str(i), ',', num2str(j));
                    obj.cdfs(key1) = trunc1;
                    
                    % Compute CDF for each state for goal 2
                    g2 = obj.goals{2};
                    mu2 = atan2(g2(2)- x(2), g2(1) - x(1)); 
                    pd2 = makedist('Normal','mu',mu2,'sigma',obj.sigma2);
                    trunc2 = truncate(pd2, -pi, pi);
                    
                    key2 = strcat(num2str(2), ',', num2str(i), ',', num2str(j));
                    obj.cdfs(key2) = trunc2;
                end
            end
            
            % Enumerate all the controls
            % UP=1, RIGHT=2, DOWN=3, LEFT=4,
            % UP_RIGHT=5, DOWN_RIGHT=6, DOWN_LEFT=7, UP_LEFT=8
            obj.controls = [1,2,3,4,5,6,7,8]; 
           
        end
        
        %% Mega prediction loop for up to timestep H.
        function preds = predict(obj, x0, H)
            
            % Convert into (i,j) index
            coor = obj.realToSim(x0);
            
            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            preds = cell([1,H+1]);
            preds(:) = {zeros(obj.gridDims)};
            
            % At current timestep, the measured state has 
            % probability = 1, zeros elsewhere.
            preds{1}(coor(1), coor(2)) = 1;
            
            for t=2:H+1
                fprintf('Predicting for t=%d\n', t);
                for g = 1:length(obj.goals)
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
                                    pug * obj.prior(g) *  preds{t-1}(s(1), s(2));
                            end
                        end
                        
                    end
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
        function prob = Pu_given_x_g(obj, u, s0, goal)
            
            if goal == 1
                key1 = strcat(num2str(goal), ',', num2str(s0(1)), ',', num2str(s0(2)));
                trunc1 = obj.cdfs(key1);
                prob = obj.cdfOfU(trunc1, u);
            elseif goal == 2
                key2 = strcat(num2str(goal), ',', num2str(s0(1)), ',', num2str(s0(2)));
                trunc2 = obj.cdfs(key2);
                prob = obj.cdfOfU(trunc2, u);
            else
                error("In PuGivenGoal(): goal is invalid: %d\n", goal);
            end
            
        end
        
        %% Converts from real state (m) to coordinate (i,j)
        function  coor = realToSim(obj, x)
            coor = [round((x(2) - obj.gridMin(2))/obj.dy) + 1, ...
                    round((x(1) - obj.gridMin(1))/obj.dx) + 1];
        end
        
        %% Converts from coordinate (i,j) to real state (m)
        function x = simToReal(obj, coor)
            x = [obj.gridMin(1) + (coor(2) - 1)*obj.dx, ...
                 obj.gridMin(2) + (coor(1) - 1)*obj.dy];
        end
        
        %% Computes the cdf of the normal distribution with mean mu and ...
        %  standard deviation sigma, evaluated at the values in the interval.
        function prob = cdfOfU(obj, truncatedDist, u)
            p = [0; 0];
            if u == 1 % UP
                p = cdf(truncatedDist,[(3*pi)/8, (5*pi)/8]);
                %p = normcdf([(3*pi)/8, (5*pi)/8],mu,sigma);
            elseif u == 2 % RIGHT
                p = cdf(truncatedDist,[-pi/8, pi/8]);
                %p = normcdf([-pi/8, pi/8],mu,sigma);
            elseif u == 3 % DOWN
                p = cdf(truncatedDist,[-(5*pi)/8, -(3*pi)/8]);
                %p = normcdf([-(5*pi)/8, -(3*pi)/8],mu,sigma);
            elseif u == 4 % LEFT
                pa = cdf(truncatedDist,[(7*pi)/8, pi]);
                pb = cdf(truncatedDist,[-pi, -(7*pi)/8]);
                p = [pa(1) + pb(1), pa(2) + pb(2)];
                %p = normcdf([(7*pi)/8, pi],mu,sigma) + ...
                %    normcdf([-pi, -(7*pi)/8],mu,sigma);
            elseif u == 5 % UP RIGHT
                p = cdf(truncatedDist,[pi/8, (3*pi)/8]);
                %p = normcdf([pi/8, (3*pi)/8],mu,sigma);
            elseif u == 6 % DOWN RIGHT
                p = cdf(truncatedDist,[-(3*pi)/8, -pi/8]);
                %p = normcdf([-(3*pi)/8, -pi/8],mu,sigma);
            elseif u == 7 % DOWN LEFT
                p = cdf(truncatedDist,[-(7*pi)/8, -(5*pi)/8]);
                %p = normcdf([-(7*pi)/8, -(5*pi)/8],mu,sigma);
            elseif u == 8 % UP LEFT
                p = cdf(truncatedDist,[(5*pi)/8, (7*pi)/8]);
                %p = normcdf([(5*pi)/8, (7*pi)/8],mu,sigma);
            end
            
            prob = p(2)-p(1);
        end
        
        %% Converts from fake discrete controls into theta controls
        function theta = uToTheta(obj, u)
            if u == 1 % UP
                theta = pi/2;
            elseif u == 2 % RIGHT
                theta = 0;
            elseif u == 3 % DOWN
                theta = -pi/2;
            elseif u == 4 % LEFT
                theta = pi;
            elseif u == 5 % UP RIGHT
                theta = pi/4;
            elseif u == 6 % DOWN RIGHT
                theta = -pi/4;
            elseif u == 7 % DOWN LEFT
                theta = -(3*pi)/4;
            elseif u == 8 % UP LEFT
                theta = (3*pi)/4;
            end
        end
        
        %% Dynamics function gives next state we can get to given current state.
        function [snext, isValid] = dynamics(obj, s0, u)
            isValid = true;
            
            if u == 1 % UP
                snext = [s0(1)+1; s0(2)];
            elseif u == 2 % RIGHT
                snext = [s0(1); s0(2)+1];
            elseif u == 3 % DOWN
                snext = [s0(1)-1; s0(2)];
            elseif u == 4 % LEFT
                snext = [s0(1); s0(2)-1];
            elseif u == 5 % UP RIGHT
                snext = [s0(1)+1; s0(2)+1];
            elseif u == 6 % DOWN RIGHT
                snext = [s0(1)-1; s0(2)+1];
            elseif u == 7 % DOWN LEFT
                snext = [s0(1)-1; s0(2)-1];
            elseif u == 8 % UP LEFT
                snext = [s0(1)+1; s0(2)-1];
            end
            
            if snext(1) < 1 || snext(1) > obj.gridDims(1) || ...
                    snext(2) < 1 || snext(2) > obj.gridDims(2)  
                snext = s0;
                isValid = false;
            end
        end
        
    end
end

