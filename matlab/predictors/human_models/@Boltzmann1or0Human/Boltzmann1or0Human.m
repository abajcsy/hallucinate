classdef Boltzmann1or0Human < DynSys
    properties
    	% Fixed human speed
        v

        % Max possible heading control bounds
        uRange
        
        % Beta (0-rand or 1-Boltzmann) fixed value
        beta

        % Probability threshold for determining likely controls
        uThresh

        % HMM model parameter
        gamma

        % Goal location for feedback policy
        theta
        
        % Time intervals between inputs
        deltaT
        
        % Num discrete controls to consider
        numCtrls

        % Number of dimensions
        dims
        
        % Store the likelyCtrls and dynamics for all states.
        likelyCtrls
        likelyMasks
        xdot
        
        % Stores discretized control increment and discretized controls.
        ctrlIncr
        discCtrls
        
        % Mixing parameter and mixing distribution 
        alpha 
        DeltaB0
        
        % (string) are we using static or dynamic beta model?
        betaModel
    end
    
    methods
        function obj = Boltzmann1or0Human(x, v, beta, uRange, gamma, ...
                theta, delta_t, uThresh, numCtrls, betaModel, extraArgs)
          %% obj = GaussianHuman(x, v, uRange, gamma, K, m, sigma, uThresh, betaModel)
          %     Dynamics of the GaussianHuman
          %         \dot{x}_1 = v * cos(u)
          %         \dot{x}_2 = v * sin(u)
          %         \dot{x}_3 = \dot{P}_t(beta = 1)
          %         -uRange(1) <= u <= uRange(2)
          %
          %     State space of the GaussianHuman
          %         x_1 = p_x
          %         x_2 = p_y
          %         x_3 = P(beta = 1)
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
          
          obj.beta = beta;

          obj.theta = theta;
          obj.deltaT = delta_t;

          obj.v = v;
          obj.gamma = gamma;
          
          obj.uRange = uRange;
          obj.uThresh = uThresh;
          obj.numCtrls = numCtrls;
          obj.ctrlIncr = (obj.uRange(2) - obj.uRange(1))/obj.numCtrls; 
          obj.discCtrls = zeros(1, obj.numCtrls);
          
          % Creates numCtrls discrete controls from (uMin, uMax]. 
          for i=1:obj.numCtrls
              obj.discCtrls(i) = obj.uRange(1) + i*obj.ctrlIncr;
          end
          
          obj.nx = length(x);
          obj.nu = 1;

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
              error("No support for beta model %s\n", betaModel);
          end
          
        end
        
        function intControls = sumDiscControls(obj, x)
            %% Approximate integral with summation: 
            %    \int_{U} e^{-\| (x_t + \Delat t f(x_t,u_t)) - \theta \|_2}

            intControls = 0.0;
            parfor i = 1:obj.numCtrls
                
                % Get discrete control.
                u = obj.discCtrls(i);
                
                % Compute the Q-value of each state and control.
                qval = obj.qFunction(x, u);
                
                % Calculate value in summation: 
                %   e^{-||(x_t + \Deltat t f(x_t,u_t)) - \theta||_2}
                val = exp(-1 .* qval);
                
                % Add to running value of summation
                intControls = intControls + obj.ctrlIncr*val; 
            end
        end
        
        function pb = betaPosterior(obj, x, u)
            %% Computes posterior given x and u
            %       P(beta=1 | xt=x, ut=u) \propto P(u | x, beta=1) * P(beta=1)
            %
            %  Note that our third state is x(3) = P(\beta=1 | xt-1, ut-1)
            
            % Posterior update in the valid range of beta
            numerator = x{3};
            
            % Approximation of integral of e^Q(x,u) over space of controls
            intControls = obj.sumDiscControls(x);
            
            % Compute the Q-value of the current control at the current
            % state.
            qval = obj.qFunction(x, u);
            pugivenx_beta1 = exp(-1 .* qval);
            
            % Get probability of other beta value (beta = 0).
            pbeta0 = 1-x{3};
            
            % Get probability of u given beta = 0 --> uniform: 1/numCtrls
            pugivenx_beta0 = 1. / obj.numCtrls;
            
            % Compute denominator of posterior.
            denominator = x{3} + ((pbeta0 .* intControls .* pugivenx_beta0) ./ (pugivenx_beta1)) ;
            
            pb = numerator ./ denominator;
            pb = max(min(pb, 1.), 0.);
            
            % Account for the probability outside the valid range
            pb = (pb .* (x{3} >= 0) .* (x{3} <= 1)) + ...
                        (x{3} .* (x{3} < 0)) + ...
                        (x{3} .* (x{3} > 1));
        end
        
        function [likelyCtrls, likelyMasks] = getLikelyControls(obj, x)
            %% Gets the set of controls that are more likely than uThresh
            %                   P(u_t | x_t) >= uThresh
            %  Input: 
            %       x           -- (cell arr) discretized states in each
            %                                 dimension
            %  Output: 
            %       likelyCtrls -- (cell arr) valid controls at each state    
            %       likelyMasks -- (Map) map with keys as the controls and 
            %                       values are a binary mask of size of the
            %                       statespace indicating if a control is
            %                       likely enough at a given state.
           
            likelyCtrls = cell(1, obj.numCtrls); % Contain all likely controls
            likelyMasks = containers.Map; % Map for likely control (str) to boolean matrix
            
            if obj.beta == 1
                %  Ground-truth human is BOLTZMANN.
                for i = 1:obj.numCtrls
                    % Compute the normalizer for the Boltzmann model
                    intControls = obj.sumDiscControls(x);

                    % Grab current candidate control.
                    u = obj.discCtrls(i);

                    % Get the Qvalue of the current state and control.
                    qval = obj.qFunction(x, u);
                    val = exp(-1 .* qval);

                    % Compute the probability.
                    pu = val ./ intControls;
                    pu = max(min(pu, 1.), 0.);

                    % Mask of same dimension as x{1}, which is 1 if coordinate is likely, 0 otherwise
                    likelyMasks(num2str(u)) = (pu >= obj.uThresh); 
                    likelyCtrls{i} = u; % Consider all controls as likely controls
                end
            else
                %  Ground-truth human is RANDOM.
                for i = 1:obj.numCtrls
                    u = obj.discCtrls(i);
                    
                    prand = 1./obj.numCtrls;   
                    
                    % Create a cell array the same size as x{1} and filled 
                    % with uniform probability over each control.
                    pu = (x{1} .* 0.) + prand;
                    pu = max(min(pu, 1.), 0.);

                    % Mask of same dimension as x{1}, which is 1 if coordinate is likely, 0 otherwise
                    likelyMasks(num2str(u)) = (pu >= obj.uThresh);
                    likelyCtrls{i} = u; % Consider all controls as likely controls
                end
            end
        end
        
        function computeUAndXDot(obj, x)
            %% Computes and stores the likley state-dependant control and state deriv.
            % Get the likely state-dependant control for the ith discrete control: 
            %   P(u_i | x)
            [obj.likelyCtrls, obj.likelyMasks] = obj.getLikelyControls(x);

            % Compute and store the corresponding dynamics.
            obj.xdot = {};
            for i=1:obj.numCtrls
                u = obj.likelyCtrls{i};
                currLikelyMask = obj.likelyMasks(num2str(u));
            	f = obj.dynamics(1,x,u); % note: the time (first arg) isnt used
                
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
        
        %% Computes the Q-function for the Boltzmann model.
        %       Q(x,u) = ||x + dt*f(x,u) - theta||_2
        function qval = qFunction(obj, x, u)
            % Find next x by forward euler
            x1 = x{1} + obj.deltaT .* obj.v .* cos(u);
            x2 = x{2} + obj.deltaT .* obj.v .* sin(u);

            % Evaluate distance of next x to goal theta under L2 norm
            qval = ((x1 - obj.theta(1)).^2 + (x2 - obj.theta(2)).^2).^(0.5);
        end
        
        %% Plots the state-dependent control distribution for Boltzmann.
        function plotPUGivenXBeta(obj, grid)
            figure(2);
            m = floor(obj.numCtrls/2);
            n = obj.numCtrls - m;
            for i=1:obj.numCtrls
                
                % Compute the normalizer for the Boltzmann model
                intControls = obj.sumDiscControls(grid.xs);

                % Grab current candidate control.
                u = obj.discCtrls(i);

                % Get the Qvalue of the current state and control.
                qval = obj.qFunction(grid.xs, u);
                pu = exp(-1 .* qval);

                % Compute the probability.
                pugivenx_beta1 = pu ./ intControls;
                pugivenx_beta1 = max(min(pugivenx_beta1, 1.), 0.);
                
                % Plot.
                subplot(m,n,i);
                pcolor(grid.xs{1}(:,:,1), grid.xs{2}(:,:,1), pugivenx_beta1(:,:,1));
                title_str = strcat('P(u=',num2str(u), '| x, beta=1)');
                title(title_str);
                colormap('bone')
                colorbar
            end
        end
    end
end

