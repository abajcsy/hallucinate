function hamValue = boltzmannBetaHuman_ham(t, data, deriv, schemeData)

if ~iscell(deriv)
  deriv = num2cell(deriv);
end

% Grab the GaussianHuman so we can access its parameters.
human = schemeData.dynSys;
x = schemeData.grid.xs;
uMode = schemeData.uMode;

%  Compute the Hamiltonian over all states and discretized controls.
pdot_f = 0.0;
parfor i=1:(2 + (length(human.betas) - 1))
    p_i = repmat(deriv{i}, [1, 1, 1, human.numCtrls]);
    pdot_f = pdot_f + (p_i .* human.xdot{i});
end
% p1 = repmat(deriv{1}, [1, 1, 1, human.numCtrls]);
% p2 = repmat(deriv{2}, [1, 1, 1, human.numCtrls]);
% p3 = repmat(deriv{3}, [1, 1, 1, human.numCtrls]);
% pdot_f = p1 .* human.xdot{1} + p2 .* human.xdot{2} + p3 .* human.xdot{3};

% Maximize/minimize Hamiltonian.
if strcmp(uMode, 'min')
	hamValue = min(pdot_f, [], 4);
else
	hamValue = max(pdot_f, [], 4);
end

% Negate hamValue if backward reachable set.
if strcmp(schemeData.tMode, 'backward')
   hamValue = -hamValue;
end

end
% function hamValue = boltzmannHuman_ham(t, data, deriv, schemeData)
% 
% if ~iscell(deriv)
%   deriv = num2cell(deriv);
% end
% 
% % Grab the GaussianHuman so we can access its parameters.
% human = schemeData.dynSys;
% x = schemeData.grid.xs;
% uMode = schemeData.uMode;
% 
% %  Compute the Hamiltonian over all states and discretized controls.
% p1 = repmat(deriv{1}, [1, 1, 1, human.numCtrls]);
% p2 = repmat(deriv{2}, [1, 1, 1, human.numCtrls]);
% p3 = repmat(deriv{3}, [1, 1, 1, human.numCtrls]);
% pdot_f = p1.*human.xdot{1} + p2.*human.xdot{2} + p3.*human.xdot{3};
% 
% % Maximize/minimize Hamiltonian
% if strcmp(uMode, 'min')
% 	hamValue = min(pdot_f, [], 4);
% else
% 	hamValue = max(pdot_f, [], 4);
% end
% 
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
% end

% function hamValue = boltzmannHuman_ham(t, data, deriv, schemeData)
% 
% if ~iscell(deriv)
%   deriv = num2cell(deriv);
% end
% 
% % Grab the GaussianHuman so we can access its parameters.
% human = schemeData.dynSys;
% x = schemeData.grid.xs;
% uMode = schemeData.uMode;
% 
% % %  Compute the Hamiltonian over all states and discretized controls.
% % p1 = repmat(deriv{1}, [1, 1, 1, human.numCtrls]);
% % p2 = repmat(deriv{2}, [1, 1, 1, human.numCtrls]);
% % p3 = repmat(deriv{3}, [1, 1, 1, human.numCtrls]);
% % pdot_f = p1.*human.xdot{1} + p2.*human.xdot{2} + p3.*human.xdot{3};
% % 
% % % Maximize/minimize Hamiltonian
% % if strcmp(uMode, 'min')
% % 	hamValue = min(pdot_f, [], 4);
% % else
% % 	hamValue = max(pdot_f, [], 4);
% % end
% 
% for i=1:human.numCtrls
%     % Get dynamics at each state given current control.
%     %xdot = human.xdot{i};
%     u = human.likelyCtrls{i};
%     xdot = human.dynamics(x,u);
%     
%     % Compute Hamiltonian: 
%     %   \grad_x V(x,t) * f(x,u).
%     hamU = deriv{1} .* xdot{1} + deriv{2} .* xdot{2} + deriv{3} .* xdot{3};
%     
%     % Maintain the min or max Hamiltonian value.
%     if strcmp(uMode, 'max')
%         hamValue = max(hamValue, hamU);
%     else
%         hamValue = min(hamValue, hamU);
%     end
% end
% end