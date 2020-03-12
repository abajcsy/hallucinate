function dx = dynamics(obj, t, x, u)
    % dx = dynamics(obj, x, u)
    %     Dynamics of the GaussianTwoGoalHuman
    %         \dot{x}_1 = v * cos(u)
    %         \dot{x}_2 = v * sin(u)
    %         \dot{x}_3 = (betaPosterior(goal = g1 | x, u) - betaPrior(goal = g1))*gamma    

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
    end
    
    dx = cell(obj.nx, 1);

    for i = 1:length(obj.dims)
        if i == 1
            % NOTE: we need to "freeze" the dynamics when we have invalid
            % probabilities. 
            dx{i} = obj.v .* cos(u) .* (x{3} >= 0) .* (x{3} <= 1);
        elseif i == 2
            dx{i} = obj.v .* sin(u) .* (x{3} >= 0) .* (x{3} <= 1);
        elseif i == 3
            if strcmp(obj.betaModel, 'static')
                % NOTE: These dynamics are for a STATIC beta model. 
                dx{i} = obj.gamma * (obj.betaPosterior(x, u) - x{3}) .* ...
                    (x{3} >= 0) .* (x{3} <= 1);
            else
                % NOTE: These dynamics are for a DYNAMIC beta model. 
                dx{i} = obj.gamma * ((obj.betaPosterior(x, u) - x{3}) + ...
                    (obj.alpha * x{3} + (1 - obj.alpha) * obj.DeltaB0 - x{3})) .* ...
                    (x{3} >= 0) .* (x{3} <= 1);
            end
        else
            error('Only dimension 1-3 are defined for dynamics of GaussianG1orG2Human!')    
        end
    end
    
%     % test: what if we freezed dynamics at the goal?
%     % NOTE: this wont work if target set is the goal region...
%     for i = 1:length(obj.dims)
%         statesInGoal = (x{1} - obj.goals{obj.trueGoalIdx}(1)).^2 + ...
%                     (x{2} - obj.goals{obj.trueGoalIdx}(2)).^2 <= obj.goalSetRad^2;
%         goalIndicies = find(statesInGoal > 0.);
%         dx{i}(goalIndicies) = 0.;
%     end
    
    if returnVector
      dx = cell2mat(dx);
    end
end