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
      currLikelyMask = likelyMasks(num2str(UCurrent));
      
      % wherever the mask is 0, change to NaN so we don't freeze
      % dynamics unecessarily. 
      currLikelyMask = currLikelyMask * 1; % convert to double arr.
      currLikelyMask(currLikelyMask == 0) = nan;
                
      hCurrent = deriv{1} .* obj.v .* cos(UCurrent) + ...
                    deriv{2} .* obj.v .* sin(UCurrent) + ...
                    deriv{3} .* (obj.gamma * (obj.betaPosterior(y, UCurrent) - y{3}));
                
      % pick out likely enough controls.          
      hCurrent = hCurrent * currLikelyMask;
      if hCurrent > hOpt
          hOpt = hCurrent;
          uOpt = UCurrent;
      end
    end

elseif strcmp(uMode, 'min')

    hOpt = 1e6;
    for i=1:length(likelyCtrls)
      UCurrent = likelyCtrls{i};
      currLikelyMask = likelyMasks(num2str(UCurrent));
      
      % wherever the mask is 0, change to NaN so we don't freeze
      % dynamics unecessarily. 
      currLikelyMask = currLikelyMask * 1; % convert to double arr.
      currLikelyMask(currLikelyMask == 0) = nan;
      
      hCurrent = deriv{1}.*obj.v.*cos(UCurrent) + ...
                deriv{2}.*obj.v.*sin(UCurrent) + ...
                deriv{3}.*(obj.gamma * (obj.betaPosterior(y, UCurrent) - y{3}));
            
      % pick out likely enough controls.     
      hCurrent = hCurrent * currLikelyMask;     
      if hCurrent < hOpt
          hOpt = hCurrent;
          uOpt = UCurrent;
      end
    end
else
  error('Unknown uMode!')
end