classdef OptimizedLatticeBayesPredictor
    %OPTIMIZEDLATTICEBAYESPREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Lattice properties.
        latticeMin
        latticeMax
        sideLength
        numRows
        numCols
        
        numStates
        numControls
        
        goals
        goalPriors
        truncDists
                
        T % Transition matrix.
        U % State-conditioned control probabilities (indexed by goal).
    end
    
    methods
        function obj = OptimizedLatticeBayesPredictor(latticeMin, latticeMax, ...
                sideLength, goals, sigmas, goalPriors)
            obj.latticeMin = latticeMin;
            obj.latticeMax = latticeMax;
            obj.sideLength = sideLength;
                        
            % Number of cells along y (i.e. number of rows).
            obj.numRows = round((obj.latticeMax(2) - obj.latticeMin(2)) / ...
                (0.5 * obj.sideLength * sqrt(3)));
            
            % Number of cells along x (i.e. number of columns).
            obj.numCols = round((obj.latticeMax(1) - obj.latticeMin(1)) / ...
                obj.sideLength);            
            
            obj.numStates = obj.numRows * obj.numCols;
            
            obj.numControls = 6;
            
            obj.goals = goals;
            
            obj.truncDists = {};
            for goalIdx = 1:length(obj.goals)
                pd = makedist('Normal', 'mu', 0, 'sigma', sigmas(goalIdx));
                obj.truncDists{goalIdx} = truncate(pd, -pi, pi);
            end
            
            obj.goalPriors = goalPriors;
            
            fprintf('Computing transition matrix...');
            obj.T = obj.computeTransitionMatrix();
            fprintf('done.\n');
            obj.U = {};
            for goalIdx = 1:length(obj.goals)
                fprintf('Computing control probability matrix for goal %d...', ...
                        goalIdx);
                obj.U{goalIdx} = obj.computeControlProbMatrix(goalIdx);
                fprintf('done.\n');
            end
        end
        
        function preds = predict(obj, x0, y0, horizon)
            [i0, j0] = obj.realToSim(x0, y0);
            s0 = (i0 - 1) * obj.numCols + j0;
            
            predsGoal = {};
            for goalIdx = 1:length(obj.goals)
                predsGoal{goalIdx} = zeros(horizon, obj.numStates, 'single');
                predsGoal{goalIdx}(1, s0) = 1;
            end
            
            parfor goalIdx = 1:length(obj.goals)
                for t = 1:(horizon - 1)
                    predsGoal{goalIdx}(t+1, :) = predsGoal{goalIdx}(t, :) * ...
                        squeeze(sum(reshape(repelem(obj.U{goalIdx}, 1, obj.numStates), ...
                        [obj.numStates, obj.numStates, obj.numControls]) .* obj.T, 3));
                end
            end
                        
            preds = zeros(horizon, obj.numStates);
            for goalIdx = 1:length(obj.goals)
                preds = preds + predsGoal{goalIdx} * obj.goalPriors(goalIdx);
            end
        end
        
        function U = computeControlProbMatrix(obj, goalIdx)
            % Compute the state-conditioned control probability matrix.
            U = zeros(obj.numStates, obj.numControls, 'single');
            
            for i = 1:obj.numRows
                for j = 1:obj.numCols
                    for u = 1:obj.numControls
                        % Fill in the state-conditioned control matrix
                        % entry.
                        s = (i - 1) * obj.numCols + j;
                        U(s, u) = obj.controlProb(u, i, j, goalIdx);
                    end
                end
            end
        end
        
        function T = computeTransitionMatrix(obj)
            % Compute the transition matrix.
            T = zeros(obj.numStates, obj.numStates, obj.numControls, 'single');
            
            for i = 1:obj.numRows
                for j = 1:obj.numCols
                    for u = 1:obj.numControls                        
                        % Fill in the transition matrix entry.
                        [inext, jnext] = obj.dynamics(i, j, u);
                        if inext >= 1 && inext <= obj.numRows && ...
                                jnext >= 1 && jnext <= obj.numCols
                            s = (i - 1) * obj.numCols + j;
                            snext = (inext - 1) * obj.numCols + jnext;
                            T(s, snext, u) = 1;
                        end
                    end
                end
            end
        end
        
        function [inext, jnext] = dynamics(obj, i, j, u)
            if u == 1 % UP_RIGHT
                if mod(i, 2) == 0
                    inext = i + 1;
                    jnext = j;
                else
                    inext = i + 1;
                    jnext = j + 1;
                end
            elseif u == 2 % RIGHT
                inext = i;
                jnext = j + 1;
            elseif u == 3 % DOWN_RIGHT
                if mod(i, 2) == 0
                    inext = i - 1;
                    jnext = j;
                else
                    inext = i - 1;
                    jnext = j + 1;
                end
            elseif u == 4 % DOWN_LEFT
                if mod(i, 2) == 0
                    inext = i - 1;
                    jnext = j - 1;
                else
                    inext = i - 1;
                    jnext = j;
                end
            elseif u == 5 % LEFT
                inext = i;
                jnext = j - 1;
            elseif u == 6 % UP_LEFT
                if mod(i, 2) == 0
                    inext = i + 1;
                    jnext = j - 1;
                else
                    inext = i + 1;
                    jnext = j;
                end
            end
        end
        
        function prob = controlProb(obj, u, i, j, goalIdx)
            % Get the lower and upper bounds to integrate the Gaussian 
            % PDF over.
            bounds = obj.getIntegrationBounds(u);
            [x, y] = obj.simToReal(i, j);
            
            g = obj.goals{goalIdx};
            mu = atan2(g(2) - y, g(1) - x);
            
            bound1 = wrapToPi(mu - bounds(1));
            bound2 = wrapToPi(mu - bounds(2));
            
            p = cdf(obj.truncDists{goalIdx}, [bound1, bound2]);
            prob = abs(p(2) - p(1));
            
            if abs(bound1 - bound2) > pi
                prob = 1 - prob;
            end
        end
        
        function bounds = getIntegrationBounds(obj, u)
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
        
        %% Converts from real state (x, y) to coordinate (i,j)
        function [i, j] = realToSim(obj, x, y)
            i = round((2 * (y - obj.latticeMin(2))) / ...
                (obj.sideLength * sqrt(3))) + 1;
            
            if mod(i, 2) == 1
                % Case where i is odd.
                j = round((x - obj.latticeMin(1)) / obj.sideLength) + 1;
            else
                % Case where i is even.
                j = round((x - obj.latticeMin(1) - 0.5 * obj.sideLength) ...
                    / obj.sideLength) + 2;
            end
        end
        
        %% Converts from coordinate (i,j) to real state (x, y)
        function [x, y] = simToReal(obj, i, j)            
            y = obj.latticeMin(2) + ...
                (i - 1) * (obj.sideLength * 0.5 * sqrt(3));
            
            if mod(i, 2) == 1
                x = obj.latticeMin(1) + (j - 1) * obj.sideLength;
            else
                x = obj.latticeMin(1) + 0.5 * obj.sideLength + ...
                    obj.sideLength * (j - 2);
            end
        end
        
        function plotPredsAtTime(obj, preds, t, eps, plotThresh)
            xs = [];
            ys = [];
            ps = [];    

            for i = 1:obj.numRows
                for j = 1:obj.numCols
                    [x, y] = obj.simToReal(i, j);
                    xs = [xs, x];
                    ys = [ys, y];
                    
                    s = (i - 1) * obj.numCols + j;

                    if ~plotThresh
                        ps = [ps, preds(t, s)];                        
                    else
                        if preds(t, s) > eps
                            ps = [ps, 1];
                        else
                            ps = [ps, 0];
                        end
                    end
                end
            end

            sum(ps)

            figure(1);
            scatter(xs, ys, 30 * ones(1, length(ys)), 1 - ps, 'filled', 'MarkerEdgeColor', 'k');
            colormap('gray');            
        end
    end
end

