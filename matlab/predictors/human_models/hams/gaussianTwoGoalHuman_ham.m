function hamValue = gaussianTwoGoalHuman_ham(t, data, deriv, schemeData)

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

% Grab the GaussianHuman so we can access its parameters.
human = schemeData.dynSys;
x = schemeData.grid.xs;
uMode = schemeData.uMode;

%  Compute the Hamiltonian over all states and discretized controls.
p1 = repmat(deriv{1}, [1, 1, 1, human.numCtrls]);
p2 = repmat(deriv{2}, [1, 1, 1, human.numCtrls]);
p3 = repmat(deriv{3}, [1, 1, 1, human.numCtrls]);
pdot_f = p1.*human.xdot{1} + p2.*human.xdot{2} + p3.*human.xdot{3};

% Maximize/minimize Hamiltonian
if strcmp(uMode, 'min')
	hamValue = min(pdot_f, [], 4);
else
	hamValue = max(pdot_f, [], 4);
end

end