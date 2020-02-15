function dx = dynamics(obj, x, u)
    % dx = dynamics(obj, x, u)
    %     Dynamics of the GaussianHuman
    %         \dot{x}_1 = v * cos(u)
    %         \dot{x}_2 = v * sin(u)
    %         \dot{x}_3 = (betaPosterior(beta=0 | x, u) - betaPrior(beta=0))*gamma    

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
    end
    
    dx = cell(obj.nx, 1);
    
    % NOTE: All dx are multiplied elementwise by likelyMasks at the
    % respective control, so as to not account for unlikely controls

    for i = 1:length(obj.dims)
        if ~isKey(obj.likelyMasks,num2str(u))
            disp("")
        end
        if i == 1
            % NOTE: we need to "freeze" the dynamics when we have invalid
            % probabilities. 
            dx{i} = obj.v .* cos(u) .* (x{3} >= 0) .* (x{3} <= 1) .* obj.likelyMasks(num2str(u));
        elseif i == 2
            dx{i} = obj.v .* sin(u) .* (x{3} >= 0) .* (x{3} <= 1) .* obj.likelyMasks(num2str(u));
        elseif i == 3
            if strcmp(obj.betaModel, 'static')
                % NOTE: These dynamics are for a STATIONARY beta model. 
                dx{i} = obj.gamma * (obj.betaPosterior(x, u) - x{3}) .* (x{3} >= 0) .* (x{3} <= 1) .* obj.likelyMasks(num2str(u));
            else
                % NOTE: These dynamics are for a DYNAMIC beta model. 
                dx{i} = obj.gamma * ((obj.betaPosterior(x, u) - x{3}) + ...
                    (obj.alpha * x{3} + (1 - obj.alpha) * obj.DeltaB0 - x{3})).* (x{3} >= 0) .* (x{3} <= 1);
            end
        else
            error('Only dimension 1-3 are defined for dynamics of GaussianHuman!')    
        end
    end

    if returnVector
      dx = cell2mat(dx);
    end
end