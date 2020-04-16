classdef GaussianTwoGoalHuman < DynSys
    %GAUSSIANTWOGOALHUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Fixed human speed
        v

        % Max possible heading control bounds
        uRange

        % Probability threshold for determining likely controls
        uThresh

        % HMM model parameter
        gamma

        % Cell array with XY locations of two (known) goals
        goals
        
        % Mean and variance in normal distribution (same for goal1 and goal2)
        sigma
        
        % Truncated Gaussian with zero mean and sigma std dev
        truncpd
        
        % Num discrete controls to consider
        numCtrls
        
        % (arr) Discrete controls \in uRange
        us
        
        % (cell) Integration bounds for each discrete control.
        ubounds

        % Number of dimensions
        dims
        
        % Store the likelyCtrls and dynamics for all states.
        likelyCtrls
        likelyMasks
        xdot
        
        % Store the optimal controls for each of the goals
        uOptG1
        uOptG2
        
        % Mixing parameter and mixing distribution 
        alpha 
        DeltaB0
        
        % (string) are we using static or dynamic beta model?
        betaModel
        
        grid
    end
    
    methods
        function obj = GaussianTwoGoalHuman(x, v, uRange, gamma, ...
                goals, sigma, uThresh, numCtrls, betaModel, extraArgs)
         %% obj = GaussianTwoGoalHuman(x, v, uRange, gamma, goals, sigma, ...
         %          uThresh, numCtrls, betaModel)
          %     Dynamics of the GaussianTwoGoalHuman
          %         \dot{x}_1 = v * cos(u)
          %         \dot{x}_2 = v * sin(u)
          %         \dot{x}_3 = \dot{P}_t(goal = g1)
          %         -uRange(1) <= u <= uRange(2)
          %
          %     State space of the GaussianTwoGoalHuman
          %         x_1 = p_x
          %         x_2 = p_y
          %         x_3 = P(goal = g1)
          %
          % TODO: replace all the 2*pi's by uRange(2) - uRange(1)
          
          if numel(x) ~= 3
            error('Initial state does not have right dimension!');
          end

          if ~iscolumn(x)
            x = x';
          end

          obj.x = x;
          obj.xhist = obj.x;
          obj.dims = 1:3;

          obj.goals = goals;
          obj.sigma = sigma;

          obj.v = v;
          obj.uRange = uRange;
          obj.gamma = gamma;
          obj.uThresh = uThresh;
          obj.numCtrls = numCtrls;
          
          obj.nx = length(x);
          obj.nu = 1;
          
          obj.xdot = [];
          obj.uOptG1 = [];
          obj.uOptG2 = [];
          obj.grid = [];
          
          % Generate numCtrls equally spaced controls in uRange.
          obj.us = obj.genControls();
          obj.ubounds = obj.genUbounds(obj.us);
          
          % truncated gaussian with zero mean.
          pd = makedist('Normal','mu',0,'sigma',obj.sigma);
          obj.truncpd = truncate(pd, obj.uRange(1), obj.uRange(2));     

          obj.betaModel = betaModel;
          if strcmp(betaModel, 'dynamic')
              if isfield(extraArgs, 'DeltaB0')
                  obj.DeltaB0 = extraArgs.DeltaB0;
              else
                  obj.DeltaB0 = 0.5;
                  warning('Setting DeltaB0 to default: 0.5\n');
              end
              if isfield(extraArgs, 'alpha')
                  obj.alpha = extraArgs.alpha;
              else
                  obj.alpha = 0.1;
                  warning('Setting alpha to default: 0.1\n');
              end
          elseif ~strcmp(betaModel, 'static')
              error('No support for beta model %s\n', betaModel);
          end
        end
        
        function pb = betaPosterior(obj, x, u)
            %% Computes posterior given x and u
            %       P(goal=g1 | xt=x, ut=u) \propto P(u | x, goal=g1) * P(goal=g1)
            %
            %  Note that our third state is x(3) = P(goal=g1 | xt-1, ut-1)
            
            % Get probability of u under each goal.
            PuGivenG1 = obj.PuGivenGoal(u, x, 1);
            PuGivenG2 = obj.PuGivenGoal(u, x, 2);
            
            % Get the probabilities of each goal.
            PG1 = x{3};
            PG2 = 1 - x{3};
            
            numerator = PuGivenG1 .* PG1;
            denominator = (PuGivenG1 .* PG1) + (PuGivenG2 .* PG2);
            
            pb = numerator./denominator;
            pb = max(min(pb, 1), 0);
            
            % Account for the probability outside the valid range
            pb = (pb .* (x{3} >= 0) .* (x{3} <= 1)) + (x{3} .* (x{3} < 0)) + (x{3} .* (x{3} > 1));
        end
                
        %% Computes P(u | x, goal) where goal=g1 or g2
        %  The observation models are:
        %   
        %       P(u | x, goal = g1) = N(atan2(g1_y - y, g1_x - x), sigma^2)
        % 
        %  and 
        % 
        %       P(u | x, goal = g2) = N(atan2(g2_y - y, g2_x - x), sigma^2)
        % 
        %  where N stands for the normal distribution. 
        function pu = PuGivenGoal(obj, u, x, goalIdx)
            %if isempty(obj.uOptG1) || isempty(obj.uOptG2)
            %    error('Optimal controls for the observation model have not been precomputed!\n');
            %    error('Make sure to first run computeUOptGoals(x).\n');
            %end
           
            % Get the optimal control and min angular distance from current
            % control to the optimal control.
            uopt = atan2(obj.goals{goalIdx}(2)- x{2}, obj.goals{goalIdx}(1) - x{1});
            diff = abs(angdiff(u, uopt));
            
            % corner case.
            if length(obj.us) == 1
                pu = 1;
                return
            end
            
            % Note because of numerics: 
            % sometimes we get controls that are just close enough to zero 
            % but are negative.
            zero_tol = -1e-7; 
            % find indicies of all controls in the positive [0, pi] range.
            pos_idxs = find(obj.us >= zero_tol);
%             % Need to make sure we grab bounds around zero.
%             for i=1:length(obj.us)
%                 bound = obj.ubounds{i};
%                 if bound(2) >= 0 && bound(1) <= 0
%                     pos_idxs(end+1) = i;
%                 end
%             end
            %note: need to include -pi in positive controls since its = pi
            %      -pi is always first in the controls.
            pos_idxs(end+1) = 1; 

            % store probabilities.
            pu = zeros(size(u));
            
             % find the control bounds.
            for i=pos_idxs
                bounds = obj.ubounds{i};
                low_bound = bounds(1);
                up_bound = bounds(2);
                
                if low_bound < 0
                    % case around zero.
                    case_around_0_idxs = find(diff <= up_bound);
                    p = cdf(obj.truncpd, [0, up_bound]);
                    pu_curr =  abs(p(2) - p(1)) * 2;
                    pu(case_around_0_idxs) = pu_curr;
                elseif up_bound < 0
                    % case around pi/-pi.
                    case_around_pi_idxs = find(diff >= low_bound);
                    p = cdf(obj.truncpd, [low_bound, pi]);
                    pu_curr =  abs(p(2) - p(1)) * 2;
                    pu(case_around_pi_idxs) = pu_curr;
                else
                    % normal case.
                    normal_bounds_idxs = find(diff <= up_bound & diff >= low_bound);
                    p = cdf(obj.truncpd, [low_bound, up_bound]);
                    pu_curr =  abs(p(2) - p(1));
                    pu(normal_bounds_idxs) = pu_curr;
                end
            end
           
        end
        
        function [likelyCtrls, likelyMasks] = getLikelyControls(obj, x)
            %% Gets the set of controls that are more likely than uThresh
            %                   P(u_t | x_t) >= uThresh
            %  Input: 
            %       x           -- (cell arr) discretized states in each
            %                                 dimension
            %  Output: 
            %       likelyCtrls -- (cell arr) valid controls at each state    
            
            likelyCtrls = cell(1, obj.numCtrls); % Contain all likely controls
            likelyMasks = containers.Map;        % Map for likely control (str) to boolean matrix
            
            for i=1:obj.numCtrls
                % Get the current discretized ontrol.
                u = obj.us(i);
                
                % Make array of this control the size of our statespace. 
                uarr = u * ones(size(x{1}));
                
                
                % Get probability of u under each goal.
                PuGivenG1 = obj.PuGivenGoal(uarr, x, 1);
                PuGivenG2 = obj.PuGivenGoal(uarr, x, 2);
                
                % Get the probabilities of each goal.
                PG1 = x{3};
                PG2 = 1 - x{3};
                
                % Compute the probability of the controls 
                % marginalized over goals. 
                PuGivenX = (PuGivenG1 .* PG1) + (PuGivenG2 .* PG2);
                
                % Pick out the controls at the states where P >= epsilon
                likelyMasks(num2str(u)) = (PuGivenX >= obj.uThresh);
                
%                 % ------ DEBUGGING ------ %
%                 tmp = (PuGivenX >= obj.uThresh) * 1;
%                 figure(3);
%                 pcolor(obj.grid.xs{1}(:,:,1), obj.grid.xs{2}(:,:,1), tmp(:,:,10));
%                 title('Opt control to get to goal1');
%                 colormap('bone')
%                 colorbar
%                 % ------ DEBUGGING ------ %
                
                likelyCtrls{i} = uarr;
            end

        end
                
        %% Gets the control bounds for integration. 
        % Returns cell array containing a vector of lower and upper
        % integration bounds for each u in us.
        function ubounds = genUbounds(obj, us)
            ubounds = cell(1,obj.numCtrls);
            incr = (obj.uRange(2) - obj.uRange(1))/obj.numCtrls;
            for i=1:obj.numCtrls
                u = us(i);
                ubounds{i} = [wrapToPi(u - incr/2), wrapToPi(u + incr/2)]; 
            end
        end
        
        %% Generate discrete num_ctrls ranging from 
        %   u \in [urange(1), urange(2)) 
        % (note: not including upper urange bound since its -pi to pi
        % and pi and -pi are the same). 
        % Returns an array of discrete controls. 
        function us = genControls(obj)
            incr = (obj.uRange(2) - obj.uRange(1))/obj.numCtrls;
            us = zeros(1,obj.numCtrls);
            u = obj.uRange(1);
            for i=1:obj.numCtrls
                us(i) = u;
                u = u + incr;
            end
        end
        
        %% Computes and stores the likley state-dependant control and state deriv.
        % Get the likely state-dependant control for the ith discrete control: 
        %   P(u_i | x)
        function computeUAndXDot(obj, x)
            [obj.likelyCtrls, obj.likelyMasks] = obj.getLikelyControls(x);

            % Compute and store the corresponding dynamics.
            obj.xdot = {};
            for i=1:obj.numCtrls
                u = obj.us(i);
                currLikelyMask = obj.likelyMasks(num2str(u));
                uarr = obj.likelyCtrls{i};
                
                % wherever the mask is 0, change to NaN so we don't freeze
                % dynamics unecessarily. 
                currLikelyMask = currLikelyMask * 1; % convert to double arr.
                currLikelyMask(currLikelyMask == 0) = nan;
                
                % note: first arg (time) is ignored.
            	f = obj.dynamics(1,x,uarr); 
                
                % Convert into an N1 x N2 x N3 x numCtrls array
                if i == 1
                    obj.xdot = f;
                    
                    obj.xdot{1} = obj.xdot{1} .* currLikelyMask;
                    obj.xdot{2} = obj.xdot{2} .* currLikelyMask;
                    obj.xdot{3} = obj.xdot{3} .* currLikelyMask;
                else
                    obj.xdot{1} = cat(4, obj.xdot{1}, f{1} .* currLikelyMask);
                    obj.xdot{2} = cat(4, obj.xdot{2}, f{2} .* currLikelyMask);
                    obj.xdot{3} = cat(4, obj.xdot{3}, f{3} .* currLikelyMask);
                end
            end
        end        
        
        
        %% Compute the optimal control over the entire statespace for each goal.
        %   These are the means of the two Gaussian models.
        function computeUOptGoals(obj, x)
            g1 = obj.goals{1};
            g2 = obj.goals{2};
            
            % Store u*goal1(x_t) and u*goal2(x_t) for all x_t
            obj.uOptG1 = atan2(g1(2)- x{2}, g1(1) - x{1}); 
            obj.uOptG2 = atan2(g2(2)- x{2}, g2(1) - x{1}); 
            
%             % ---- Debugging.  ---- %
%             figure(3);
%             pcolor(obj.grid.xs{1}(:,:,1), obj.grid.xs{2}(:,:,1), obj.uOptG1(:,:,1));
%             title('Opt control to get to goal1');
%             colormap('bone')
%             colorbar
% 
%             figure(4);
%             pcolor(obj.grid.xs{1}(:,:,1), obj.grid.xs{2}(:,:,1), obj.uOptG2(:,:,1));
%             title('Opt control to get to goal2');
%             colormap('bone')
%             colorbar
%             % ---- Debugging.  ---- %
        end
        
        %% Helper: plots the probability of the given action at each state.
        function plotPu(obj, u, pu, goal)
            figure(goal);
            puXY = pu(:,:,1);
            pcolor(obj.grid.xs{1}(:,:,1), obj.grid.xs{2}(:,:,1), puXY);
            title(strcat('P(u=', num2str(u), ' | g=goal', num2str(goal), ')'));
            colorbar
            caxis([0,1])
        end
        
        %% Helper: set grid for debugging through plotting.
        function setGrid(obj, grid)
            obj.grid = grid;
        end
    end
end

