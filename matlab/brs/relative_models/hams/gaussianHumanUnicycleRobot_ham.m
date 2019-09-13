function hamValue = gaussianHumanUnicycleRobot_ham(t, data, deriv, schemeData)

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

% Grab the GaussianHuman so we can access its parameters.
dyn = schemeData.dynSys;
x = schemeData.grid.xs;
vRange = dyn.vRRange;
wRange = dyn.wRRange;
uMode = schemeData.uMode;
dMode = schemeData.dMode;

%  Compute the Hamiltonian over all states and discretized controls.
p1 = repmat(deriv{1}, [1, 1, 1, 1, dyn.numCtrls]);
p2 = repmat(deriv{2}, [1, 1, 1, 1, dyn.numCtrls]);
p3 = repmat(deriv{3}, [1, 1, 1, 1, dyn.numCtrls]);
pdot_f = p1.*dyn.xdot{1} + p2.*dyn.xdot{2} + p3.*dyn.xdot{3};

% Maximize/minimize Hamiltonian for disturbnaces
if strcmp(dMode, 'min')
	hamValue = min(pdot_f, [], 5);
else
	hamValue = max(pdot_f, [], 5);
end

% Maximize/minimize Hamiltonian for control
speedDecision = deriv{1} .* cos(x{4}) + deriv{2} .* sin(x{4});
if strcmp(uMode, 'min')
  hamRobot = ((speedDecision >= 0) * vRange(1)) + ...
    ((speedDecision < 0) * vRange(2)) + ((deriv{4} >= 0) * wRange(1)) + ...
    ((deriv{4} < 0) * wRange(2));
else
	hamRobot = ((speedDecision >= 0) * vRange(2)) + ...
    ((speedDecision < 0) * vRange(1)) + ((deriv{4} >= 0) * wRange(2)) + ...
    ((deriv{4} < 0) * wRange(1));
end
hamValue = hamValue + hamRobot;

% Flip the hamiltonian for BRS
hamValue = -hamValue;

% % Compare against the Hamiltonian obtained by just maximizing the first
% % two terms
% uOpt_nobeta = atan2(deriv{2}, deriv{1});
% xdot_nobeta = human.dynamics(x, uOpt_nobeta);
% hamU_nobeta = deriv{1} .* xdot_nobeta{1} + deriv{2} .* xdot_nobeta{2} + deriv{3} .* xdot_nobeta{3};
% if strcmp(uMode, 'max')
%     condition = (human.likelyCtrls{1} - uOpt_nobeta < 0) + (human.likelyCtrls{end} - uOpt_nobeta > 0) - 1e-2;
%     
%     hamValue = max(hamValue, hamU_nobeta) .* (sign(condition) <= 0) + ...
%         hamValue .* (sign(condition) > 0);
% end

end