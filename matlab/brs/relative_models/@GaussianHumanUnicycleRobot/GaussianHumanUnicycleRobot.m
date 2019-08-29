classdef GaussianHumanUnicycleRobot < DynSys
    properties
        % Fixed human speed
        v

        % Max possible heading control bounds for the human
        uRange
        
        % Speed range for the robot
        vRRange
        
        % Angular speed range for the robot
        wRRange

        % Probability threshold for determining likely controls for the
        % humans
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
        
        % Normalizer for the gaussian distribution
        gaussianNorm
        
        % Store the likelyCtrls and dynamics for all states.
        likelyCtrls
        xdot
        
        % Mixing parameter and mixing distribution 
        alpha 
        DeltaB0
        
        % (string) are we using static or dynamic beta model?
        betaModel
        
        
    end
    
    methods
        function obj = GaussianHumanUnicycleRobot(x, v, uRange, gamma, ...
                K, m, sigma, uThresh, numCtrls, betaModel, vRRange, ...
                wRRange, extraArgs)
          
          if numel(x) ~= 4
            error('Initial state does not have right dimension!');
          end

          if ~iscolumn(x)
            x = x';
          end

          obj.x = x;
          obj.xhist = obj.x;
          obj.dims = 1:4;

          obj.K = K;
          obj.m = m;
          obj.sigma = sigma;

          obj.v = v;
          obj.uRange = uRange;
          obj.vRRange = vRRange;
          obj.wRRange = wRRange;
          obj.gamma = gamma;
          obj.uThresh = uThresh;
          obj.numCtrls = numCtrls;

          obj.nx = length(x);
          obj.nu = 2;
          obj.nd = 1;

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
            uOpt = 0;
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
            
            % Compute optimal control: u* = Kx + m
            optCtrl = 0;
            
            % Compute the inner part that we are going to take log of.
            % TODO: Safeguard against x{3} being zero.
            innerLog = (sqrt(2*pi*obj.sigma^2))*((obj.uThresh - (1/(2*pi))*(1 - x{3}))./x{3});
            
            % NOTE: We need to safeguard against cases where we are taking 
            %       a square-root of a negative number. To do this, we 
            %       need to ensure that the result of the log is > 0 and < 1. 
            C = sqrt(-2*obj.sigma^2 .* log(innerLog)) .* (innerLog > 0) .* (innerLog <= 1) + ...
                1e6 * (innerLog <= 0) + 1e6 * (innerLog > 1) + ...
                1e6 * (x{3} < 0) + 1e6 * (x{3} > 1);
            
            % Compute the bounds on the likley controls.
            upperBound = min(optCtrl + C, obj.uRange(2));
            lowerBound = max(optCtrl - C, obj.uRange(1));
            
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

            % Compute and store the corresponding dynamics corresponding to
            % the zero control of the robot. 
            obj.xdot = {};
            u = {0, 0};
            for i=1:obj.numCtrls
                d = obj.likelyCtrls{i};
                f = obj.dynamics(x, u, d);
                % Convert into an N1 x N2 x N3 x numCtrls array
                if i == 1
                    obj.xdot = {f{1}, f{2}, f{3}};
                else
                    obj.xdot{1} = cat(5, obj.xdot{1}, f{1});
                    obj.xdot{2} = cat(5, obj.xdot{2}, f{2});
                    obj.xdot{3} = cat(5, obj.xdot{3}, f{3});
                end
            end
        end
        
    end
end

