classdef GaussianHuman < DynSys
    properties
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

        % Distribution in HMM
        DeltaB0

        % Number of dimensions
        dims
    end
    
    methods
        function obj = GaussianHuman(x, v, uRange, gamma, ...
                K, m, sigma, uThresh, DeltaB0, numCtrls)
          %% obj = GaussianHuman(x, v, uRange, gamma, K, m, sigma, uThresh, DeltaB0)
          %     Dynamics of the GaussianHuman
          %         \dot{x}_1 = v * cos(u)
          %         \dot{x}_2 = v * sin(u)
          %         \dot{x}_3 = (DeltaB0(beta=0) - betaPosterior(beta=0 | x, u))*gamma
          %         -pi <= u <= pi
          %     
          %     State space of the GaussianHuman
          %         x_1 = p_x
          %         x_2 = p_y
          %         x_3 = P(beta = 0)
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

          obj.DeltaB0 = DeltaB0;
          
        end
        
        function pb = betaPosterior(obj, x, u)
            %% Computes posterior given x and u
            %       P(beta=0 | xt=x, ut=u) \propto P(u | x, beta=0) * P(beta=0)
            %
            %  Note that our third state is x(3) = P(\beta=0 | xt-1, ut-1)
            
            uOpt = obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m;
            numerator = x{3};
            denominator = x{3} + sqrt((2*pi*obj.sigma^2)/(2*pi)) .* ...
                (1 - x{3}) .* 1./exp(-(u - uOpt).^2 ./ 2*obj.sigma^2);
            
            pb = numerator./denominator;
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
            optCtrl = obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m;
            
            % Compute the inner part that we are going to take log of.
            % TODO: Safeguard against x{3} being zero.
            innerLog = (sqrt(2*pi*obj.sigma^2)./x{3}).*(obj.uThresh - (1/(2*pi)).*(1 - x{3}));
            
            % NOTE: We need to safeguard against cases where we are taking 
            %       a square-root of a negative number. To do this, we 
            %       need to ensure that the result of the log is > 0 and < 1. 
            C = sqrt(-2*obj.sigma^2 .* log(innerLog)) .* (innerLog > 0) .* (innerLog < 1) + ...
                1e6 * (innerLog < 0) + 1e6 * (innerLog > 1);
            
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
        
    end
end

