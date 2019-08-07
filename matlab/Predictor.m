classdef Predictor
    %PREDICTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xinit       % (array) start position (in grid cells)
        g           % (array) known goal location (in grid cells)
        betas       % (array) discrete values beta can take
        Pbeta       % (map) prior over beta -- beta values are keys, probability is value
        gridDims    % (array) num grid cells in each dim (i.e. height and width)
        states      % (cell arr) indicies of all states in grid
        controls    % (cell arr) all controls
    end
    
    methods
        function obj = Predictor(xinit, goal, gridDims, beta_values, beta_prior)
            %PREDICTOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.xinit = xinit;
            obj.g = goal;
            obj.betas = beta_values;
            obj.Pbeta = beta_prior; 
            obj.gridDims = gridDims;
            obj.Pbeta = containers.Map(obj.betas, beta_prior);
            
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
        function preds = predict(obj, x0, H)

            % Initialize empty prediction grids forward in time.
            % Assume P(xt | x0) = 0 for all xt
            preds = cell([1,H+1]);
            preds(:) = {zeros(obj.gridDims)};
            
            % Current measured state has probability = 1, zeros elsewhere.
            preds{1}(x0(1), x0(2)) = 1;
            
            for t=2:H+1
                for xcurr = obj.states
                    xt = xcurr{1};
                    for b = obj.betas
                        Pb = obj.Pbeta(b);
                        for ut1 = obj.controls
                            % Invert dynamics to get the state x_t-1 we came
                            % from applying control u_t-1.
                            [xt1, isValid] = obj.invDyn(xt, ut1);
                            
                            if isValid
                                % Get the probability of the prior state, x_t-1
                                Pxt1_x0 = preds{t-1}(xt1(1), xt1(2));

                                % Get the probability of this action being
                                % taken from x_t-1, given this value of beta.
                                Put1_xt1_beta = obj.Pu_given_x_b(ut1, xt1, b);

                                % Compute the probability of this new state:
                                % p(xt | x0) += P(beta)*P(ut-1|xt-1,beta)*P(xt-1|x0)
                                preds{t}(xt(1), xt(2)) = preds{t}(xt(1), xt(2)) + ...
                                                         Pb * Put1_xt1_beta * Pxt1_x0;
                            end
                        end
                    end
                end
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
            offset = 0; %max(cell2mat(values(Qs)));

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
        
        %% Plotting. 
        function plot(obj, preds, time)
            % Plot the posterior over beta.
            figure(1)
            heatmap(cell2mat(values(obj.Pbeta)), 'ColorLimits',[0 1]);
            ax = gca;
            ax.XData = cell2mat(keys(obj.Pbeta)); 
            ax.YData = ["P(b)"];

            % Plot state distributions.
            figure(2)
            %hold on
            for t = time
                heatmap(preds{t}, 'ColorLimits',[0 1], 'Colormap', cool, 'FontSize', 4);
                %im = imagesc(preds{t});
                %im.AlphaData = .6;
                %colormap(cool);
            end
            %colorbar
            %ax = gca;
            %ax.Title = "P(xH|x0)";
        end
    end
end

