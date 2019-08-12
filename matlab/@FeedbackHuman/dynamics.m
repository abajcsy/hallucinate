function dx = dynamics(obj, ~, x, u, ~)
    % dx = dynamics(obj, ~, x, u)
    %     Dynamics of the FeedbackHuman
    %         \dot{x}_1 = v * cos(u)
    %         \dot{x}_2 = v * sin(u)
    %         \dot{x}_3 = alpha * betaPosterior(b=0 | x, u) + (1-alpha) * betaPrior(b=0)

    dx = cell(obj.nx, 1);

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
    end

    for i = 1:length(obj.dims)
        if i == 1
            dx{i} = obj.v .* cos(u);
        elseif i == 2
            dx{i} = obj.v .* sin(u);
        elseif i == 3
            dx{i} = obj.alpha * obj.betaPosterior(x, u) + (1-obj.alpha) * obj.betaPrior;
        else
            error('Only dimension 1-3 are defined for dynamics of FeedbackHuman!')    
        end
    end

    if returnVector
      dx = cell2mat(dx);
    end
end
