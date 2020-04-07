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
[likelyCtrls, likelyMasks] = obj.getLikelyControls(y);

if strcmp(uMode, 'max')
    hOpt = -1e6;
    for i=1:length(likelyCtrls)
      UCurrent = likelyCtrls{i};
      currentMask = likelyMasks(num2str(UCurrent));
      currentMask = currentMask * 1;
      currentMask(currentMask == 0) = nan;
      hCurrent = (deriv{1}.*obj.v.*cos(UCurrent) + ...
                deriv{2}.*obj.v.*sin(UCurrent) + ...
                deriv{3}.*(obj.gamma * (obj.betaPosterior(y, UCurrent) - y{3}))).*currentMask;
      if hCurrent > hOpt
          hOpt = hCurrent;
          uOpt = UCurrent;
      end
    end

elseif strcmp(uMode, 'min')

    hOpt = 1e6;
    for i=1:length(likelyCtrls)
      UCurrent = likelyCtrls{i};
      currentMask = likelyMasks(num2str(UCurrent));
      currentMask = currentMask * 1;
      currentMask(currentMask == 0) = nan;
      hCurrent = (deriv{1}.*obj.v.*cos(UCurrent) + ...
                deriv{2}.*obj.v.*sin(UCurrent) + ...
                deriv{3}.*(obj.gamma * (obj.betaPosterior(y, UCurrent) - y{3}))).*currentMask;
      if hCurrent < hOpt
          hOpt = hCurrent;
          uOpt = UCurrent;
      end
    end
else
  error('Unknown uMode!')
end