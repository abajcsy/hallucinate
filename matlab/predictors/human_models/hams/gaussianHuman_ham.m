function hamValue = gaussianHuman_ham(t, data, deriv, schemeData)

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

% Grab the GaussianHuman so we can access its parameters.
human = schemeData.dynSys;
x = schemeData.grid.xs;
uMode = schemeData.uMode;

% Get the current range of likely controls. 
likelyCtrls = human.getLikelyControls(x);

if strcmp(uMode, 'min')
  hamValue = Inf;
elseif strcmp(uMode, 'max')
  hamValue = -Inf; 
else
 error('Unknown uMode!')
end

for i=1:human.numCtrls
    % Get the likely state-dependant control for the ith discrete control: 
    %   P(u_i | x)
    u = likelyCtrls{i};
    
    % Get dynamics at each state given current control.
    xdot = human.dynamics(x,u);
    
    % Compute Hamiltonian: 
    %   \grad_x V(x,t) * f(x,u).
    hamU = deriv{1} .* xdot{1} + deriv{2} .* xdot{2} + deriv{3} .* xdot{3};
    
    % Maintain the min or max Hamiltonian value.
    if strcmp(uMode, 'max')
        hamValue = max(hamValue, hamU);
    else
        hamValue = min(hamValue, hamU);
    end
end

end