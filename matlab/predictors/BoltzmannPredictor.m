classdef BoltzmannPredictor
    %PREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xinit       % (array) start position (in grid cells)
        g           % (array) known goal location (in grid cells)
        
        betas       % (array) discrete values beta can take
        priorBeta   % (map) prior over beta -- beta values are keys, probability is value
        beliefBeta  % (map) belief over beta, updated after receiving measurements
        
        gridDims    % (array) num grid cells in each dim (i.e. height and width)
        states      % (cell arr) indicies of all states in grid
        controls    % (cell arr) all controls
    end
    
    methods
        function obj = BoltzmannPredictor(xinit, goal, gridDims, beta_values, beta_prior)
            %PREDICTOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.xinit = xinit;
            obj.g = goal;
            obj.betas = beta_values;
            obj.gridDims = gridDims;
            obj.priorBeta = containers.Map(obj.betas, beta_prior);
            obj.beliefBeta = containers.Map(obj.betas, beta_prior);
            
            % enumerate all the state indicies
            obj.states = {};
            for x=1:obj.gridDims(1)
                for y=1:obj.gridDims(2)
                    obj.states{end+1} = [x;y];
                end
            end
            
            % enumerate all the controls
            % UP=1, RIGHT=2, DOWN=3, LEFT=4,
            % UP_RIGHT=5, DOWN_RIGHT=6, DOWN_LEFT=7, UP_LEFT=8
            obj.controls = [1,2,3,4,5,6,7,8]; 
        end
        
        
        %% Predicts the forward reachable sets H steps into the future given x0 
        %  Given:
        %       x0          -- initial state
        %       H           -- prediction horizon
        %       ctrlThresh  -- threshold for which controls are likely enough 
        %  Output:
        %       preds       -- cell array indexed by 1:H+1 with corresponding FRS
        function preds = predictDeterministic(obj, x0, H, ctrlThresh)

            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            preds = cell([1,H+1]);
            preds(:) = {zeros(obj.gridDims)};
            
            % At current timestep, the measured state has 
            % probability = 1, zeros elsewhere.
            preds{1}(x0(1), x0(2)) = 1;
            
            % Store belief as we "hallucinate" observations during prediction. 
            hallucinateBelief = obj.beliefBeta;
 
            statesToConsider = {x0};
            for t=2:H+1
                nextStates = {};
                for idx=1:length(statesToConsider)
                    % (0) Get the current state to compute FRS for. 
                    xt = statesToConsider{idx};
                    
                    % (1) Compute likely controls at this state.
                    controlSet = obj.getLikelyControls(xt, ctrlThresh, hallucinateBelief);

                    % (2) Compute the likeliood of getting to state xt 
                    %     given the likely controls
                    frs = obj.computeFRS(xt, controlSet);
                    for state = frs
                        xnext = state{1};
                        preds{t}(xnext(1), xnext(2)) = 1;
                    end

                    % (3) Find the worst-case control under the current
                    %     distribution over beta, P(beta)
                    worstU = obj.getWorstCaseControl(hallucinateBelief);

                    % (4) Update the hallucinated belief
                    hallucinateBelief = obj.updateBeliefBeta(worstU, xt, hallucinateBelief);

                    % Record the next states to compute FRS for. 
                    nextStates = horzcat(nextStates, frs);
                end
                statesToConsider = nextStates;
            end
        end
        
        
        %% Predicts the state distribution H steps into the future given x0 
        %  Computes:
        %       P(xt | x0) = \sum_xt-1 P(x_t | x_t-1)*P(x_t-1 | x_0)
        %                  
        %  Given:
        %       x0    -- initial state
        %       H     -- prediction horizon
        %  Output:
        %       preds -- cell array indexed by 1:H+1 with corresponding state
        %                distributions.
        function preds = predictStochastic(obj, x0, H)

            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            preds = cell([1,H+1]);
            preds(:) = {zeros(obj.gridDims)};
            
            % At current timestep, the measured state has 
            % probability = 1, zeros elsewhere.
            preds{1}(x0(1), x0(2)) = 1;
            
            % Store belief as we "hallucinate" observations during prediction. 
            hallucinateBelief = obj.betaBelief;
            
            % Threshold for which controls are likely enough to be
            % considered.
            ctrlThresh = 0;
            
            for t=2:H+1
                for xcurr = obj.states
                    xt = xcurr{1};
                    
                    % (1) compute likely controls at this state.
                    controlSet = obj.getLikelyControls(xt, ctrlThresh, hallucinateBelief);
                    
                    % (2) compute the likeliood of getting to state xt 
                    %     given the likely controls
                    preds{t}(xt(1), xt(2)) = ...
                        obj.computeStateLikelihood(xt, controlSet, hallucinateBelief, preds{t-1});
                    
                    % (3) find the worst-case control under the current
                    %     distribution over beta, P(beta)
                    worstU = obj.getWorstCaseControl(hallucinateBelief);
                    
                    % (4) update the hallucinated belief
                    hallucinateBelief = obj.updateBeliefBeta(worstU, xt, hallucinateBelief);
                end
            end
        end
        
        %% Gets the worst-case control under the curr distribution of beta.
        function worstU = getWorstCaseControl(obj, beliefBeta)
            %TODO: this is hard-coded based on my toy example!
            %      change this to do something smarter, like KL div or
            %      something. 
            worstU = 3;
        end
        
        %% Returns map of likely controls and their probabilities
        function controlSet = getLikelyControls(obj, xt, ctrlThresh, hallucinateBelief)
            controlSet = containers.Map('KeyType','int32', 'ValueType','any');
            for u = obj.controls
                prob = 0;
                for b = obj.betas
                    prob = prob + obj.Pu_given_x_b(u, xt, b)*hallucinateBelief(b);
                end
                if prob > ctrlThresh
                    controlSet(u) = prob;
                end
            end
        end
        
        %% Computes the likelihood of reaching x_t 
        %  given the  the set of possible controls in controlSet 
        %  and the current distribution over beta parameters. 
        function prob = computeStateLikelihood(obj, xt, controlSet, beliefBeta, prevStateDist)
            prob = 0;
            for b = obj.betas
                Pb = beliefBeta(b);

                for u = keys(controlSet)
                    ut1 = u{1};
                    
                    % Invert dynamics to get the state x_t-1 we came
                    % from applying control u_t-1.
                    [xt1, isValid] = obj.invDyn(xt, ut1);

                    if isValid
                        % Get the probability of the prior state, x_t-1
                        Pxt1_x0 = prevStateDist(xt1(1), xt1(2));

                        % Get the probability of this action being
                        % taken from x_t-1, given this value of beta.
                        Put1_xt1_beta = obj.Pu_given_x_b(ut1, xt1, b);

                        % Compute the probability of this new state:
                        % p(xt | x0) += P(beta)*P(ut-1|xt-1,beta)*P(xt-1|x0)
                        prob = prob + Pb * Put1_xt1_beta * Pxt1_x0;
                    end
                end
            end
        end
        
        %% Updates prior over beta via HMM from RSS paper. 
        %     P^k_-(beta) = (1-epsilon)P^k-1_+(beta) + epsilon*P^0_-(beta)
        %  Inputs:
        %       P_k1    -- (map) contains old posterior
        %       epsilon -- (float) epsilon probability that beta is
        %                   resampled from the original prior
        function P_k = P_beta_HMM(obj, P_k1, epsilon)
            P_k = P_k1;
            
            for b = obj.betas
                P_k(b) = (1-epsilon)*P_k1(b) + epsilon*obj.priorBeta(b);
            end
        end
        
        %% Compute the probability of a specific action given a state and beta value.
        %  Boltzmann observation model:
        %       P(u | x0; beta) \propto e^{-beta*Q(x0, u)}
        %  where the Q-function is simply:
        %       Q(x0, u) = ||u||_2 + ||x0 + u - g||_2
        %  and we assume that each control action is norm 1. 
        function prob = Pu_given_x_b(obj, u, x0, beta)

            % Compute the next states that we could possibly get to given our
            % dynamics. 
            Qs = containers.Map; % keys are controls, values are Q-values
            for ctrl = obj.controls
                [xnext, isValid] = obj.dynamics(x0,ctrl);
                if isValid
                    if ctrl >= 1 && ctrl <= 4
                        % straight actions have 1 cost.
                        Qx0u = 1 + norm(xnext - obj.g);
                    else
                        % diagonal actions have sqrt(2) action cost.
                        Qx0u = sqrt(2) + norm(xnext - obj.g);
                    end
                    % store Q value in map.
                    Qs(num2str(ctrl)) = Qx0u;
                end
            end

            % Normalization trick to improve numerical stability.
            % See: http://cs231n.github.io/linear-classify/#softmax
            offset = max(cell2mat(values(Qs)));

            % Compute the denominator by summing over all possible actions.
            normalizer = 0;
            for q = values(Qs)
                normalizer = normalizer + exp(-log(beta) * q{1} - offset);
            end

            % Compute the numerator by evaluating the likelihood of the given action. 
            prob = exp(-log(beta) * Qs(num2str(u)) - offset)/normalizer;
        end
        
        %% Inverts dynamics to find which state we were at previously.
        %  Given deterministic dynamics:
        %           xnext = f(xprev, u)
        %  this function solves for:
        %           xprev = finv(xnext, u)
        function [xprev, isValid]= invDyn(obj, xnext, u)
            isValid = true;
            
            if u == 1 % UP
                xprev = [xnext(1)+1; xnext(2)];
            elseif u == 2 % RIGHT
                xprev = [xnext(1); xnext(2)-1];
            elseif u == 3 % DOWN
                xprev = [xnext(1)-1; xnext(2)];
            elseif u == 4 % LEFT
                xprev = [xnext(1); xnext(2)+1];
            elseif u == 5 % UP RIGHT
                xprev = [xnext(1)+1; xnext(2)-1];
            elseif u == 6 % DOWN RIGHT
                xprev = [xnext(1)-1; xnext(2)-1];
            elseif u == 7 % DOWN LEFT
                xprev = [xnext(1)-1; xnext(2)+1];
            elseif u == 8 % UP LEFT
                xprev = [xnext(1)+1; xnext(2)+1];
            end
                
            if xprev(1) < 1 || xprev(1) > obj.gridDims(2) || ...
                    xprev(2) < 1 || xprev(2) > obj.gridDims(1)
                % NOTE: gridDims is X,Y but the coordinates in xprev = [Y, X]!!
                %       this matters when gridDims isnt the same for X Y -- fix
                %       this!
                xprev = xnext;
                isValid = false;
            end
        end
        
        %% Dynamics function gives next state we can get to given current state.
        function [xnext, isValid] = dynamics(obj, x0, u)
            isValid = true;
            
            if u == 1 % UP
                xnext = [x0(1)-1; x0(2)];
            elseif u == 2 % RIGHT
                xnext = [x0(1); x0(2)+1];
            elseif u == 3 % DOWN
                xnext = [x0(1)+1; x0(2)];
            elseif u == 4 % LEFT
                xnext = [x0(1); x0(2)-1];
            elseif u == 5 % UP RIGHT
                xnext = [x0(1)-1; x0(2)+1];
            elseif u == 6 % DOWN RIGHT
                xnext = [x0(1)+1; x0(2)+1];
            elseif u == 7 % DOWN LEFT
                xnext = [x0(1)+1; x0(2)-1];
            elseif u == 8 % UP LEFT
                xnext = [x0(1)-1; x0(2)-1];
            end
            
            if xnext(1) < 1 || xnext(1) > obj.gridDims(2) || ...
                    xnext(2) < 1 || xnext(2) > obj.gridDims(1)  
                % NOTE: gridDims is X,Y but the coordinates in xprev = [Y, X]!!
                %       this matters when gridDims isnt the same for X Y -- fix
                %       this!
                xnext = x0;
                isValid = false;
            end
        end
        
        %% Computes the set of states that are reachable in one timestep
        %  from state x_t given the control set. 
        function frs = computeFRS(obj, xt, controlSet)
            frs = {};
            for u = keys(controlSet)
                ut = u{1};
                [xnext, isValid] = obj.dynamics(xt, ut);
                if isValid
                    frs{end+1} = xnext;
                end
            end
        end
        
        %% Update posterior over beta given measured u0
        %  Apply Bayesian update:  
        %       P'(beta | x0, u0) \propto P(u0 | x0; beta)*P(beta)
        %  Given observation (u0) and the state from which this action was
        %  taken (x0).
        function priorNext = updateBeliefBeta(obj, u0, x0, beliefBeta)
            posterior = beliefBeta;
                
%             % Compute normalizer. 
%             denominator = 0;
%             for b = obj.betas
%                 denominator = denominator + ...
%                     obj.Pu_given_x_b(u0, x0, b)*beliefBeta(b);
%             end
% 
%             % Update the distribution over beta.
%             for b = obj.betas
%                 numerator = obj.Pu_given_x_b(u0, x0, b)*beliefBeta(b);
%                 posterior(b) = numerator/denominator;
%             end
            
            % Apply "epsilon-static" transition model to beta.
            epsilon = 0.02;
            priorNext = obj.P_beta_HMM(posterior, epsilon);
        end
        
        %% Plotting. 
        function plot(obj, preds, time)
            % Plot the posterior over beta.
            figure(1)
            heatmap(cell2mat(values(obj.priorBeta)), 'ColorLimits',[0 1]);
            ax = gca;
            ax.XData = cell2mat(keys(obj.priorBeta)); 
            ax.YData = ["P(b)"];

            % Plot state distributions.
            figure(2)
            for t = time
                subplot(length(time),1,t);
                heatmap(preds{t}, 'ColorLimits',[0 1], 'Colormap', cool, 'FontSize', 4);
                title(strcat('t=',num2str(t)));
                colorbar('off')
            end
        end
    end
end

