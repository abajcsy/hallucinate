classdef FeedbackHuman < DynSys
  properties
    % Fixed human speed
    v
    
    % Max possible heading control bounds
    uRange
    
    % Probability threshold for determining likely controls
    uThresh
    
    % HMM model parameter
    alpha
    
    % Gains for the feedback control policy
    K
    m
    
    % Prior over beta=0 
    betaPrior
    
    % Number of dimensions
    dims
  end
  
  methods
    function obj = FeedbackHuman(x, v, uRange, alpha, K, m, uThresh, prior)
      %% obj = FeedbackHuman(x, v, uRange, alpha)
      %     Dynamics of the FeedbackHuman
      %         \dot{x}_1 = v * cos(u)
      %         \dot{x}_2 = v * sin(u)
      %         \dot{x}_3 = (betaPrior(beta=0) - betaPosterior(beta=0 | x, u))*alpha
      %         -pi <= u <= pi
      %     
      %     State space of the FeedbackHuman
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
      
      obj.v = v;
      obj.uRange = uRange;
      obj.alpha = alpha;
      obj.uThresh = uThresh;
      
      obj.nx = length(x);
      obj.nu = 1;
      
      obj.betaPrior = prior;
    end
    
    function pb = betaPosterior(obj, x, u)
        %% Computes posterior given x and u
        %       P(beta=0 | xt=x, ut=u) \propto P(u | x, beta=0) * P(beta=0)
        %
        %   Note that since P(u | x, beta=0) = 1, and our third state
        %       x(3) = P(\beta=0 | xt-1, ut-1)
        %   then we have:
        %       P(beta=0 | xt=x, ut=u) \propto P(\beta=0 | xt-1, ut-1) \propto x(3)
        
        if ~iscell(x)
          x = num2cell(x);
        end
        
        if u == obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m
            pb = x{3}./obj.P_u_given_x(u, x);
        else
            pb = 0;
        end
    end
    
    function pux = P_u_given_x(obj, u, x)
        %% Computes the distribution over actions given x and 
        %   marginalized over the betas
        %       P(u | x) = \sum_beta P(u | x, beta) * P(beta)
        
        if ~iscell(x)
            x =  num2cell(x);
        end
        
        if u == obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m
            pux = x{3} + (1 - x{3})/2*pi;
        else
            pux = (1 - x{3})/2*pi;
        end
    end
    
    function validURange = getLikelyControls(obj, x)
        %% Gets the set of controls that are more likely than uThresh
        if (1/2*pi)*(1-obj.betaPrior) > obj.uThresh
            validURange = obj.uRange;
        else
            validURange = obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m;
        end
    end
    
    
  end % end methods
end % end classdef
