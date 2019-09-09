classdef Beta1or0Human < DynSys
    %BETA1OR0HUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Beta (0 or 1) fixed value
        beta
        
        % Fixed human speed
        v

        % Max possible heading control bounds
        uRange

        % Probability threshold for determining likely controls
        uThresh

        % HMM model parameter
        gamma
        
        % Gains for the feedback control policy
        K
        m
        
        % Variance in normal distribution
        sigma
        
        % Num discrete controls to consider
        numCtrls

        % Number of dimensions
        dims
        
        % Store the likelyCtrls and dynamics for all states.
        likelyCtrls
        xdot
        
        % Mixing parameter and mixing distribution 
        alpha 
        DeltaB0
        
        % (string) are we using static or dynamic beta model?
        betaModel
        
        % Normalizer for the gaussian distribution
        gaussianNorm
    end
    
    methods
        function obj = Beta1or0Human(x, v, beta, uRange, gamma, ...
                K, m, sigma, uThresh, numCtrls, betaModel, extraArgs)
          %% obj = Beta1or0Human(x, v, beta, uRange, gamma, uThresh,numCtrls, betaModel)
          %     Dynamics of the GaussianHuman
          %         \dot{x}_1 = v * cos(u)
          %         \dot{x}_2 = v * sin(u)
          %         \dot{x}_3 = \dot{P}_t(beta = 0)
          %         -uRange(1) <= u <= uRange(2)
          %
          %     State space of the RandomHuman
          %         x_1 = p_x
          %         x_2 = p_y
          %         x_3 = P(beta = 0)
          %
          %     Human policy model
          %         beta = 0 --> u ~ N(Kx+m, sigma^2)
          %         beta = 1 --> u ~ rand(uRange(1), uRange(2))
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
          
          obj.K = K;
          obj.m = m;
          obj.sigma = sigma;

          obj.v = v;
          obj.uRange = uRange;
          obj.gamma = gamma;
          obj.uThresh = uThresh;
          obj.numCtrls = numCtrls;

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
         
          mu = obj.K*obj.x(1:2) + obj.m;
          if mu ~= 0
              error('This code only works for K=0 and m=0!');
          end
          
          % Normalize the gaussian to take into account control bounds.
          p = normcdf([obj.uRange(1), obj.uRange(2)], mu, obj.sigma);
          obj.gaussianNorm = p(2)-p(1);
        end
        
        function pb = betaPosterior(obj, x, u)
            %% Computes posterior given x and u
            %       P(beta=0 | xt=x, ut=u) \propto P(u | x, beta=0) * P(beta=0)
            %
            %  Note that our third state is x(3) = P(\beta=0 | xt-1, ut-1)
            
            % Posterior update in the valid range of beta
            uOpt = obj.K(1)*x{1} + obj.K(2)*x{2} + obj.m;
            numerator = x{3};

            denominator = x{3} + ((sqrt(2*pi*obj.sigma^2)*obj.gaussianNorm)/(2*pi)) * ...
                (1 - x{3}) .* (exp(((u - uOpt).^2) / (2*obj.sigma^2)));
            
            pb = numerator./denominator;
            pb = max(min(pb, 1), 0);
            
            % Account for the probability outside the valid range
            pb = (pb .* (x{3} >= 0) .* (x{3} <= 1)) + (x{3} .* (x{3} < 0)) + (x{3} .* (x{3} > 1));
        end
        
        function likelyCtrls = getLikelyControls(obj, x)
            %% Gets the set of controls that are more likely than uThresh
            %                   P(u_t | x_t) >= uThresh
            %  Input: 
            %       x           -- (cell arr) discretized states in each
            %                                 dimension
            %  Output: 
            %       likelyCtrls -- (cell arr) valid controls at each state    
            
            %TBD
            if obj.beta == 0 % u drawn from Normal Distribution
                upperBound = 0*x{3};
                lowerBound = 0*x{3};
                
                % Compute optimal control: u* = Kx + m
                optCtrl = obj.K(1)*x{1} + obj.K(2)*x{2} + obj.m;
            
                % Compute the inner part that we are going to take log of.
                innerLog = sqrt(2*pi*obj.sigma^2)*obj.uThresh;
                
                % NOTE: We need to safeguard against cases where we are taking 
                %       a square-root of a negative number. To do this, we 
                %       need to ensure that the result of the log is > 0 and < 1. 
                C = sqrt(-2*obj.sigma^2 .* log(innerLog)) .* (innerLog > 0) .* (innerLog <= 1) + ...
                    1e6 * (innerLog <= 0) + 1e6 * (innerLog > 1);
            
                % TODO: need to make this size of state space???
                
                % Compute the bounds on the likley controls.
                upperBound = min(upperBound + optCtrl + C, obj.uRange(2));
                lowerBound = max(lowerBound + optCtrl - C, obj.uRange(1));
            else % u drawn from Rand(umin, umax)
                upperBound = 0*x{3};
                lowerBound = 0*x{3};
                if (1/(2*pi)) >= obj.uThresh 
                    lowerBound = obj.uRange(1)+lowerBound;
                    upperBound = obj.uRange(2)+upperBound;
                end
            end
            
            % Based on N number of discrete contrls, partition controls
            % between lower and upper bound state-wise
            likelyCtrls = cell(1, obj.numCtrls);
            linNums = linspace(0,1,obj.numCtrls);
            parfor i=1:obj.numCtrls
                likelyCtrls{i} = linNums(i)*lowerBound + (1-linNums(i))*upperBound;
            end
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
            	f = obj.dynamics(x,u);
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
    end
end

