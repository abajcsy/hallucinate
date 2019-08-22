function dx = dynamics(obj, x, u)
    % dx = dynamics(obj, x, u)
    %     Dynamics of the GaussianHuman
    %         \dot{x}_1 = v * cos(u)
    %         \dot{x}_2 = v * sin(u)
    %         \dot{x}_3 = (betaPrior(beta=0) - betaPosterior(beta=0 | x, u))*alpha    

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
    end
    
    dx = cell(obj.nx, 1);

    for i = 1:length(obj.dims)
        if i == 1
            dx{i} = obj.v .* cos(u);
        elseif i == 2
            dx{i} = obj.v .* sin(u);
        elseif i == 3
            %dx{i} = (obj.DeltaB0 - obj.betaPosterior(x, u)) * obj.gamma;
            
            % NOTE: These dynamics are for a STATIONARY beta model. 
            dx{i} = (obj.betaPosterior(x, u) - x{3}) * obj.gamma;
        else
            error('Only dimension 1-3 are defined for dynamics of GaussianHuman!')    
        end
    end

    if returnVector
      dx = cell2mat(dx);
    end
end