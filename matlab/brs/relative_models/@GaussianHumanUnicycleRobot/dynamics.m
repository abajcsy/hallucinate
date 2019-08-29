function dx = dynamics(obj, x, u, d)
    % dx = dynamics(obj, x, u)
    %     Dynamics of the GaussianHuman and Unicycle Robot.
    %         \dot{x}_1 = v_R * cos(x{3}) - v_H * cos(u_H)
    %         \dot{x}_2 = v_R * sin(x{3}) - v_H * sin(u_H)
    %         \dot{x}_3 = (betaPosterior(beta=0 | x, u) - betaPrior(beta=0))*gamma
    %         \dot{x}_4 = omega_R

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
    end
    
    if ~iscell(u)
      u = num2cell(u);
    end
    
    dx = cell(obj.nx, 1);

    for i = 1:length(obj.dims)
        if i == 1
            % NOTE: we need to "freeze" the dynamics when we have invalid
            % probabilities. 
            dx{i} = (u{1} .* cos(x{4}) - obj.v * cos(d)) .* (x{3} >= 0) .* (x{3} <= 1);
        elseif i == 2
            dx{i} = (u{1} .* sin(x{4}) - obj.v * sin(d)) .* (x{3} >= 0) .* (x{3} <= 1);
        elseif i == 3
            if strcmp(obj.betaModel, 'static')
                % NOTE: These dynamics are for a STATIONARY beta model. 
                dx{i} = obj.gamma * (obj.betaPosterior(x, d) - x{3}) .* (x{3} >= 0) .* (x{3} <= 1);
            else
                % NOTE: These dynamics are for a DYNAMIC beta model. 
                dx{i} = obj.gamma * ((obj.betaPosterior(x, d) - x{3}) + ...
                    (obj.alpha * x{3} + (1 - obj.alpha) * obj.DeltaB0 - x{3})).* (x{3} >= 0) .* (x{3} <= 1);
            end
        elseif i == 4
            dx{i} = u{2} .* (x{3} >= 0) .* (x{3} <= 1);
        else
            error('Only dimension 1-3 are defined for dynamics of GaussianHuman!')    
        end
    end

    if returnVector
      dx = cell2mat(dx);
    end
end