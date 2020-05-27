classdef BoltzmannBetaHuman < DynSys
    properties
    	% Fixed human speed
        v

        % Max possible heading control bounds
        uRange

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
        
        % List of beta parameters
        betas
    end
    
    methods
        function obj = BoltzmannBetaHuman(x, v, uRange, gamma, ...
                betas, theta, delta_t, uThresh, numCtrls, betaModel, extraArgs)
          %% obj = GaussianHuman(x, v, uRange, gamma, K, m, sigma, uThresh, betaModel)
          %     Dynamics of the GaussianHuman
          %         \dot{x}_1 = v * cos(u)
          %         \dot{x}_2 = v * sin(u)
          %         \dot{x}_3 = \dot{P}_t(beta = 0)
          %         -uRange(1) <= u <= uRange(2)
          %
          %     State space of the GaussianHuman
          %         x_1 = p_x
          %         x_2 = p_y
          %         x_3 = P(beta = 0)
          %
          % TODO: replace all the 2*pi's by uRange(2) - uRange(1)
          
          if numel(x) ~= (2 + (length(betas) - 1))
            error('Initial state does not have right dimension!');
          end

          if ~iscolumn(x)
            x = x';
          end

          obj.x = x;
          obj.xhist = obj.x;
          obj.dims = 1:(2 + (length(betas) - 1));
          
          obj.theta = theta;
          obj.betas = betas;
          obj.deltaT = delta_t;

          obj.v = v;
          obj.uRange = uRange;
          obj.gamma = gamma;
          obj.uThresh = uThresh;
          obj.numCtrls = numCtrls;
          obj.ctrlIncr = (obj.uRange(2) - obj.uRange(1))/obj.numCtrls; 
          obj.discCtrls = zeros(1, obj.numCtrls);
          for i = 1:obj.numCtrls
              obj.discCtrls(i) = obj.uRange(1) + obj.ctrlIncr*i;
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
          
          % Normalize the gaussian to take into account control bounds.
%           p = normcdf([obj.uRange(1), obj.uRange(2)], mu, obj.sigma);
%           obj.gaussianNorm = p(2)-p(1);
          
        end
        
%         function intControls = sumDiscControls(obj, x)
%             %% Approximate integral with summation: 
%             %    \int_{U} e^{-\| (x_t + \Delat t f(x_t,u_t)) - \theta \|_2}
%             % Note maybe use PUgivenBeta function?
% 
% 
%         end
        
%         function pb = betaPosterior(obj, x, u, beta_index)
%             %% Computes posterior given x and u
%             %       P(beta=1 | xt=x, ut=u) \propto P(u | x, beta=1) * P(beta=1)
%             %
%             %  Note that our third state is x(3) = P(\beta=1 | xt-1, ut-1)
%             
%             beta = obj.betas{beta_index};
%             prior = x{2+beta_index};
%             
%             pb = (obj.PUGivenBeta(u,x,beta) .* prior) ./ obj.PUGivenX(u,x);
%             pb = max(min(pb, 1.), 0.);
%             
%             % Account for the probability outside the valid range
%             pb = (pb .* (x{2+beta_index} >= 0) .* (x{2+beta_index} <= 1)) + ...
%                         (x{2+beta_index} .* (x{2+beta_index} < 0)) + ...
%                         (x{2+beta_index} .* (x{2+beta_index} > 1));
%         end
        
        function betaPosteriors = betaPosterior(obj, x, u)
            %% Computes posterior given x and u
            %       P(beta=1 | xt=x, ut=u) \propto P(u | x, beta=1) * P(beta=1)
            %
            %  Note that our third state is x(3) = P(\beta=1 | xt-1, ut-1)
            betaPosteriors = cell(1,length(obj.betas));
            sumPosteriors = 0.0;
            for i=1:length(obj.betas)-1
                beta = obj.betas{i};
                prior = x{2+i};

                pb = (obj.PUGivenBeta(u,x,beta) .* prior) ./ obj.PUGivenX(u,x);
                pb = max(min(pb, 1.), 0.);

                % Account for the probability outside the valid range
                betaPosteriors{i} = (pb .* (x{2+i} >= 0) .* (x{2+i} <= 1)) + ...
                        (x{2+i} .* (x{2+i} < 0)) + ...
                        (x{2+i} .* (x{2+i} > 1));
                    
                sumPosteriors = sumPosteriors + betaPosteriors{i};
            end
            
            % Make sure that sum of probabilites are inside the valid range
            for i=1:length(obj.betas)-1
                betaPosteriors{i} = (betaPosteriors{i} .* (sumPosteriors >= 0) .* (sumPosteriors <= 1)) + ...
                        (x{2+i} .* (sumPosteriors < 0)) + ...
                        (x{2+i} .* (sumPosteriors > 1));
            end
        end
        
        function pb = PUGivenBeta(obj, u, x, beta)
            intControls = 0.0;
            numerator = 0.0;
            parfor i = 1:obj.numCtrls
                
                % Get discrete control.
                u_i = obj.discCtrls(i);
                
                % Compute the Q-value of each state and control.
                qval = obj.qFunction(x, u_i);
                
                % Calculate value in summation: 
                %   e^{-||(x_t + \Deltat t f(x_t,u_t)) - \theta||_2}
                val = exp(beta .* qval);
                
                % Add to running value of summation - note not
                % approximating intergral so no ctrlIncr*val.
                intControls = intControls + val; 
            end
            numerator = exp(beta .* obj.qFunction(x, u));
            pb = numerator ./ intControls;
        end
        
        function pb = PUGivenX(obj,u,x)
            pb = 0.0;
            prior_total = 0.0; % running total of priors to calculate prior of last beta (which is not explicitly in the state space
            for j=1:length(obj.betas)
                beta = obj.betas{j};
                if j < length(obj.betas)
                    prior = x{j+2};
                    prior_total = prior_total + prior;
                else
                    prior = 1 - prior_total; % find efficient way of doing 1 - x{3} - ...
                end
                pb = pb + (obj.PUGivenBeta(u, x, beta) .* prior);
            end
            pb = max(min(pb, 1), 0);
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
            likelyMasks = containers.Map; % Map for likely control (str) to boolean matrix
            for i = 1:obj.numCtrls
                % Grab current candidate control.
                u = obj.discCtrls(i);
                
                % Find probability of u
                pb = obj.PUGivenX(u,x);
                likelyMasks(num2str(u)) = pb >= obj.uThresh; % Mask of same dimension as x{1}, which is 1 if coordinate is likely, 0 otherwise
                likelyCtrls{i} = u; % Consider all controls as likely controls
            end
        end
        
        %% Computes the Q-function for the Boltzmann model.
        %       Q(x,u) = ||x + dt*f(x,u) - theta||_2
        function qval = qFunction(obj, x, u)
            % Find next x by forward euler
            x1 = x{1} + obj.deltaT .* obj.v .* cos(u);
            x2 = x{2} + obj.deltaT .* obj.v .* sin(u);

            % Evaluate distance of next x to goal theta under L2 norm
            qval = -((x1 - obj.theta(1)).^2 + (x2 - obj.theta(2)).^2).^(0.5);
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
                
                currLikelyMask = currLikelyMask * 1;
                currLikelyMask(currLikelyMask == 0) = nan;
                
                % Convert into an N1 x N2 x N3 x numCtrls array
                if i == 1
                    obj.xdot = f;
                    for j=1:(2 + (length(obj.betas) - 1))
                        obj.xdot{j} = obj.xdot{j} .* currLikelyMask;
                    end
                else
                    for j=1:(2 + (length(obj.betas) - 1))
                        obj.xdot{j} = cat((2 + length(obj.betas)), obj.xdot{j}, f{j} .* currLikelyMask);
                    end
                end
            end
        end
        
    end
end

