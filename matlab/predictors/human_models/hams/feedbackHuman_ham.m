function hamValue = feedbackHuman_ham(t, data, deriv, schemeData)

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

%% Optimal control
uMode = schemeData.uMode;
K = schemeData.dynSys.K;
m = schemeData.dynSys.m;
v = schemeData.dynSys.v;
alpha = schemeData.dynSys.alpha;
obj = schemeData.dynSys;

% Get the current range of likely controls. 
validURange = obj.getLikelyControls(schemeData.grid.xs);

% Get the optimal theta for each state. 
thetaUOpt = K(1)*schemeData.grid.xs{1} + K(2)*schemeData.grid.xs{2} + m; 

% Compute the hamiltonian associated with the optimal theta. 
hOptU = deriv{1}.*v.*cos(thetaUOpt) + deriv{2}.*v.*sin(thetaUOpt) + ...
       deriv{3}.*((obj.DeltaB0 - obj.betaPosterior(schemeData.grid.xs, thetaUOpt))*alpha);

% To test with prior, not posterior
% hOptU = deriv{1}.*v.*cos(thetaUOpt) + deriv{2}.*v.*sin(thetaUOpt) + ...
%         deriv{3}.*((obj.DeltaB0 - schemeData.grid.xs{3})*alpha);
     
if strcmp(uMode, 'max')

    hOptMax = sqrt(deriv{1}.*deriv{1} + deriv{2}.*deriv{2})*v + ...
        deriv{3}.*(obj.DeltaB0*alpha);

% To test with prior, not posterior
%     hOptMax = sqrt(deriv{1}.*deriv{1} + deriv{2}.*deriv{2})*v + ...
%         deriv{3}.*((obj.DeltaB0 - schemeData.grid.xs{3})*alpha);
    
    hamValue = hOptU .* (validURange{2} - validURange{1} <= 0) + ...
        (validURange{2} - validURange{1} > 0) .* ...
        (hOptMax .* (hOptMax - hOptU >= 0) + hOptU .* (hOptMax - hOptU < 0));

elseif strcmp(uMode, 'min')

    % lambda1 * v * cos(theta) + lambda2 * v * sin(theta)
    % is minimized at atan(lambda2*v / lambda1*v) - pi
    hOptMin = -sqrt(deriv{1}.*deriv{1} + deriv{2}.*deriv{2})*v + ...
        deriv{3}.*(obj.DeltaB0*alpha);

    hamValue = hOptMin .* (hOptMin - hOptU <= 0) .* (validURange{2} - validURange{1} > 0) + ...
        hOptU .* (hOptMin - hOptU > 0).* (validURange{2} - validURange{1} <= 0);

else
  error('Unknown uMode!')
end






end