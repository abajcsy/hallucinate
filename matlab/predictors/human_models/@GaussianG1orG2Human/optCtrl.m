function uOpt = optCtrl(obj, t, y, deriv, uMode)
% uOpt = optCtrl(obj, t, y, deriv, uMode)

%% Input processing
if nargin < 5
  uMode = 'min';
end

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

if ~iscell(y)
  y = num2cell(y);
end

%% Optimal control

% Get the current range of likely controls. 
likelyCtrls = obj.getLikelyControls(y);

% % Get the optimal control for each state for GROUND TRUTH GOAL.
% goal = obj.goals{obj.trueGoalIdx};
% thetaUOpt = atan2(goal(2)- y{2}, goal(1) - y{1});
% 
% % see if the set of likely controls contains the optimal one.
% %index = find([likelyCtrls{:}] == thetaUOpt);
% 
% % if ~isempty(index)
% %     uOpt = thetaUOpt;
% % else
% 
% % Compute the hamiltonian associated with the optimal theta. 
% hOptU = deriv{1} .* obj.v .* cos(thetaUOpt) + ...
%         deriv{2} .* obj.v .* sin(thetaUOpt) + ...
%         deriv{3} .* (obj.gamma * (obj.betaPosterior(y, thetaUOpt) - y{3}));
% 
% lambda1 = deriv{1}.*obj.v;
% lambda2 = deriv{2}.*obj.v;
% thetaStar = atan2(lambda2,lambda1);
% 
% if strcmp(uMode, 'max')
% 
%     % lambda1 * v * cos(theta) + lambda2 * v * sin(theta)
%     % is maximized at atan(lambda2*v / lambda1*v)
%     hOptMax = deriv{1}.*obj.v.*cos(thetaStar) + ...
%               deriv{2}.*obj.v.*sin(thetaStar) + ...
%               deriv{3}.*(obj.gamma * (obj.betaPosterior(y, thetaStar) - y{3}));
% 
%     if hOptMax > hOptU
%         uOpt = thetaStar;
%     else
%         uOpt = thetaUOpt;
%     end
% 
% elseif strcmp(uMode, 'min')
% 
%     % lambda1 * v * cos(theta) + lambda2 * v * sin(theta)
%     % is minimized at atan(lambda2*v / lambda1*v) - pi
%     hOptMin = deriv{1}.*obj.v.*cos(thetaStar- pi) + ...
%               deriv{2}.*obj.v.*sin(thetaStar- pi) + ...
%               deriv{3}.*(obj.gamma * (obj.betaPosterior(y, thetaStar) - y{3}));
% 
%     if hOptMin < hOptU
%         uOpt = wrapToPi(thetaStar - pi);
%     else
%         uOpt = thetaUOpt;
%     end
% 
% else
%   error('Unknown uMode!')
% end

if strcmp(uMode, 'max')
    hOpt = -1e6;
    for i=1:length(likelyCtrls)
      UCurrent = likelyCtrls{i};
      hCurrent = deriv{1}.*obj.v.*cos(UCurrent) + ...
                deriv{2}.*obj.v.*sin(UCurrent) + ...
                deriv{3}.*(obj.gamma * (obj.betaPosterior(y, UCurrent) - y{3}));
      if hCurrent > hOpt
          hOpt = hCurrent;
          uOpt = UCurrent;
      end
    end

elseif strcmp(uMode, 'min')

    hOpt = 1e6;
    for i=1:length(likelyCtrls)
      UCurrent = likelyCtrls{i};
      hCurrent = deriv{1}.*obj.v.*cos(UCurrent) + ...
                deriv{2}.*obj.v.*sin(UCurrent) + ...
                deriv{3}.*(obj.gamma * (obj.betaPosterior(y, UCurrent) - y{3}));
      if hCurrent < hOpt
          hOpt = hCurrent;
          uOpt = UCurrent;
      end
    end
else
  error('Unknown uMode!')
end