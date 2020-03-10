%% Grab the sequence of optimal controls starting from a state x. 
%  Returns cell array of all the optimal controls.
function uopt = GetOptControls(x, grid, valueFuns, times, human, uMode)
    uopt = cell(1,length(times));
    for t=1:length(times)
        % Grab the derivative at all states.
        deriv = computeGradients(grid, valueFuns(:,:,:,t));

        % Value of the derivative at that particular state
        current_deriv = eval_u(grid, deriv, x);

        % Get the optimal control to apply at this state
        u = human.optCtrl(t, x, current_deriv, uMode); 
        uopt{t} = u;
    end
end