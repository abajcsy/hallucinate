function uOpt = optCtrl(obj, t, y, deriv, uMode)
% uOpt = optCtrl(obj, t, y, deriv, uMode)

%% Input processing
if nargin < 5
  uMode = 'min';
end

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

%% Optimal control

% Get the current range of likely controls. 
validURange = obj.getLikelyControls(y);

% Get the optimal theta for each state. 
thetaUOpt = obj.K(1)*y{1} + obj.K(2)*y{2} + obj.m; 

if iscell(validURange)
    uOpt = thetaUOpt;
else
    % Compute the hamiltonian associated with the optimal theta. 
    hOptU = deriv{1}.*obj.v.*cos(thetaUOpt) + deriv{2}.*obj.v.*sin(thetaUOpt) + ...
            deriv{3}.*(obj.alpha*obj.betaPosterior(y, thetaUOpt) + (1-obj.alpha)*obj.betaPrior);

    lambda1 = deriv{1}.*obj.v;
    lambda2 = deriv{2}.*obj.v;
    thetaStar = atan2(lambda2,lambda1);

    if strcmp(uMode, 'max')

        % lambda1 * v * cos(theta) + lambda2 * v * sin(theta)
        % is maximized at atan(lambda2*v / lambda1*v)
        hOptMax = deriv{1}.*obj.v.*obj.cos(thetaStar) + deriv{2}.*obj.v.*sin(thetaStar) + ...
            deriv{3}.*((1-obj.alpha) * obj.betaPrior);

        if hOptMax > hOptU
            uOpt = thetaStar;
        else
            uOpt = thetaUOpt;
        end

    elseif strcmp(uMode, 'min')

        % lambda1 * v * cos(theta) + lambda2 * v * sin(theta)
        % is minimized at atan(lambda2*v / lambda1*v) - pi
        hOptMin = deriv{1}.*obj.v.*cos(thetaStar- pi) + deriv{2}.*obj.v.*sin(thetaStar- pi) + ...
            deriv{3}.*((1-obj.alpha) * obj.betaPrior);

        if hOptMin < hOptU
            uOpt = wrapToPi(thetaStar - pi);
        else
            uOpt = thetaUOpt;
        end

    else
      error('Unknown uMode!')
    end

end





end