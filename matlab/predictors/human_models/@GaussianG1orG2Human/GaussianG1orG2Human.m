classdef GaussianG1orG2Human < DynSys
    %GAUSSIANG1ORG2HUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Ground-truth value of goal (can = 1 or 2)
        trueGoalIdx
        
        % Radius around goal (defining goal region)
        goalSetRad
        
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
        
        % Num discrete controls to consider
        numCtrls

        % Number of dimensions
        dims
        
        % Store the likelyCtrls and dynamics for all states.
        likelyCtrls
        xdot
        
        % Store the optimal controls for each of the goals
        uOptG1
        uOptG2
        
        % Mixing parameter and mixing distribution 
        alpha 
        DeltaB0
        
        % (string) are we using static or dynamic beta model?
        betaModel
        
        % Normalizer for the gaussian distribution
        gaussianNorm
        
        grid
    end
    
    methods
        function obj = GaussianG1orG2Human(x, v, trueGoalIdx, goalSetRad, uRange, ...
                gamma, goals, sigma, uThresh, numCtrls, betaModel, ...
                extraArgs)
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
          %     Human policy model
          %         trueGoalIdx = 1 --> u ~ N(u*(g1), sigma^2)
          %         trueGoalIdx = 2 --> u ~ N(u*(g2), sigma^2)
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
          obj.trueGoalIdx = trueGoalIdx;
          obj.goalSetRad = goalSetRad;

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
        function pu = PuGivenGoal(obj, u, x, goal)
            if isempty(obj.uOptG1) || isempty(obj.uOptG2)
                error('Optimal controls for the observation model have not been precomputed!\n');
                error('Make sure to first run computeUOptGoals(x).\n');
            end   
            c0 = 1/(sqrt(2*pi*obj.sigma^2));
            if length(u) ~= 1
                if goal == 1
                    errorG1 = wrapToPi(u - obj.uOptG1);
                    uG1Diff = -(errorG1).^2;
                    pu = c0 .* exp(uG1Diff ./ (2*obj.sigma^2));
                elseif goal == 2
                    errorG2 = wrapToPi(u - obj.uOptG2);
                    uG2Diff = -(errorG2).^2;
                    pu = c0 .* exp(uG2Diff ./ (2*obj.sigma^2));
                else
                    error('In PuGivenGoal(): goal is invalid: %d\n', goal);
                end
            else
                if goal == 1
                    g1 = obj.goals{1};
                    uOptForG1 = atan2(g1(2)- x{2}, g1(1) - x{1}); 
                    errorG1 = wrapToPi(u - uOptForG1);
                    uG1Diff = -(errorG1).^2;
                    pu = c0 .* exp(uG1Diff ./ (2*obj.sigma^2));
                elseif goal == 2
                    g2 = obj.goals{2};
                    uOptForG2 = atan2(g2(2)- x{2}, g2(1) - x{1}); 
                    errorG2 = wrapToPi(u - uOptForG2);
                    uG2Diff = -(errorG2).^2;
                    pu = c0 .* exp(uG2Diff ./ (2*obj.sigma^2));
                else
                    error('In PuGivenGoal(): goal is invalid: %d\n', goal);
                end
            end
        end
        
        function likelyCtrls = getLikelyControls(obj, x)
            %% Gets the set of controls that are more likely than uThresh
            %                   P(u_t | x_t, goal = g*) >= uThresh
            %  Input: 
            %       x           -- (cell arr) discretized states in each
            %                                 dimension
            %  Output: 
            %       likelyCtrls -- (cell arr) valid controls at each state    
            
            blend = linspace(0,1,obj.numCtrls);
            binaryMap = cell(1, obj.numCtrls);
            candidateCtrls = cell(1, obj.numCtrls);
            for i=1:obj.numCtrls
                % Get the current discretized control. 
                u = blend(i) * (obj.uRange(1)*ones(size(x{1}))) + ...
                    (1-blend(i)) * (obj.uRange(2)*ones(size(x{1})));
                
                % Sanity check. 
                if obj.trueGoalIdx ~= 1 && obj.trueGoalIdx ~= 2
                    error('TrueGoalIdx is invalid! Must be 1 or 2.')  
                end
                
                % Get probability of u under THE GROUND TRUTH GOAL.
                %                  P(u | x, g = g*)
                PuGivenG = obj.PuGivenGoal(u, x, obj.trueGoalIdx);
                
                % ---- Debugging. ---- %
                %obj.plotPu(u(1,1,1), PuGivenG, obj.trueGoalIdx);
                % ---- Debugging. ---- %
                
                % Pick out the controls at the states where P >= epsilon
                binaryMap{i} = PuGivenG >= obj.uThresh; 
                candidateCtrls{i} = u;
            end
            
            % Concatinate all the binary maps and candidate controls
            catBinaryMap = cell2mat(reshape(binaryMap,1,1,1,[]));
            catCandidateCtrls = cell2mat(reshape(candidateCtrls,1,1,1,[]));
            
            % Multiply the two to pick out the valid controls for each
            % state
            validCtrls = catBinaryMap .* catCandidateCtrls;
            
            % NOTE: This isn't proper masking! Need to do this hack.
            minValidCtrls = validCtrls;
            minValidCtrls(find(catBinaryMap == 0)) = 10000000000.0;
            maxValidCtrls = validCtrls;
            maxValidCtrls(find(catBinaryMap == 0)) = -10000000000.0;
            
            % Take the min and max over the product in the last dimension
            lowerBound = min(minValidCtrls, [], 4);
            upperBound = max(maxValidCtrls, [], 4);
            
            % Based on N number of discrete contrls, partition controls
            % between lower and upper bound state-wise
            likelyCtrls = cell(1, obj.numCtrls);
            
            % Minimum angular distance between the control bounds. 
            diff = angdiff(lowerBound, upperBound);
            
            % Direction to traverse unit circle when linearly interpolating.
            reverse_dir = sign(diff); 
            
            % increment to add to control.
            incr = abs(diff)/(obj.numCtrls-1);

            % generate controls \in [lowerBound, upperBound]
            for i=1:obj.numCtrls
                likelyCtrls{i} = ((i-1)*incr + lowerBound) .* (reverse_dir > 0) + ...
                                    wrapToPi(lowerBound - (i-1)*incr) .* (reverse_dir < 0);
            end

%             linNums = linspace(0,1,obj.numCtrls);
%             parfor i=1:obj.numCtrls
%                 %for i=1:obj.numCtrls
%                 % NOTE: may need to take care of some angle wrapping shit
%                 % here.... linear interpolation moves counterclockwise 
%                 % but when the action is going to the left, then we need to
%                 % interpolate clockwise.
%                 likelyCtrls{i} = linNums(i)*lowerBound + (1-linNums(i))*upperBound;
%             end
        end
        
        function computeUAndXDot(obj, x)
            %% Computes and stores the likley state-dependant control and state deriv.
            % Get the likely state-dependant control for the ith discrete control: 
            %   P(u_i | x)
            obj.likelyCtrls = obj.getLikelyControls(x);

            % Compute and store the corresponding dynamics.
            obj.xdot = {};
            for i=1:obj.numCtrls
                u = obj.likelyCtrls{i};
            	f = obj.dynamics(1,x,u); % note: first arg (time) is ignored.
                % Convert into an N1 x N2 x N3 x numCtrls array
                if i == 1
                    obj.xdot = f;
                else
                    obj.xdot{1} = cat(4, obj.xdot{1}, f{1});
                    obj.xdot{2} = cat(4, obj.xdot{2}, f{2});
                    obj.xdot{3} = cat(4, obj.xdot{3}, f{3});
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

