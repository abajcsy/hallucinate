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
    
    % Distribution in HMM
    DeltaB0
    
    % Number of dimensions
    dims
  end
  
  methods
    function obj = FeedbackHuman(x, v, uRange, alpha, K, m, uThresh, DeltaB0)
      %% obj = FeedbackHuman(x, v, uRange, alpha)
      %     Dynamics of the FeedbackHuman
      %         \dot{x}_1 = v * cos(u)
      %         \dot{x}_2 = v * sin(u)
      %         \dot{x}_3 = (DeltaB0(beta=0) - betaPosterior(beta=0 | x, u))*alpha
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
      
      obj.DeltaB0 = DeltaB0;
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
        
        uOpt = obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m;
        pb = (x{3}./obj.P_u_given_x(u, x)) .* (uOpt == u) + ...
            0 * (abs(uOpt - u) > 0);  
    end
    
    function pux = P_u_given_x(obj, u, x)
        %% Computes the distribution over actions given x and 
        %   marginalized over the betas
        %       P(u | x) = \sum_beta P(u | x, beta) * P(beta)
        
        if ~iscell(x)
            x =  num2cell(x);
        end
        
        uOpt = obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m;
        pux = (x{3} + (1 - x{3})/(2*pi)) .* (uOpt == u) + ...
            ((1 - x{3})/(2*pi)) .* (abs(uOpt - u) > 0);
    end
    
    function validURange = getLikelyControls(obj, x)
        %% Gets the set of controls that are more likely than uThresh
        optCtrl = obj.K(1).*x{1} + obj.K(2).*x{2} + obj.m;
        uMin = obj.uRange(1) .* ((1/(2*pi))*(1-x{3}) - obj.uThresh > 0) + ...
               optCtrl .* ((1/(2*pi))*(1-x{3}) - obj.uThresh <= 0);
        
        uMax = obj.uRange(2) .* ((1/(2*pi))*(1-x{3}) - obj.uThresh > 0) + ...
                optCtrl .* ((1/(2*pi))*(1-x{3}) - obj.uThresh <= 0);
            
        validURange = {uMin, uMax};
    end
    
    
  end % end methods
end % end classdef
