function dx = dynamics(obj, t, x, u)
    % dx = dynamics(obj, x, u)
    %     Dynamics of the GaussianHuman
    %         \dot{x}_1 = v * cos(u)
    %         \dot{x}_2 = v * sin(u)
    %         \dot{x}_3 = (betaPosterior(beta=1 | x, u) - betaPrior(beta=1))*gamma    

    returnVector = false;
    if ~iscell(x)
      returnVector = true;
      x = num2cell(x);
    end
    
    dx = cell(obj.nx, 1);
    
%     % NOTE: All dx are multiplied elementwise by likelyMasks at the
%     % respective control, so as to not account for unlikely controls
%     if isempty(obj.likelyMasks)
%         likelyMasks = 1.0; % if likelyMasks isnt specified, make everything likely. 
%     else
%         if ~isKey(obj.likelyMasks,num2str(u))
%             error(strcat(num2str(u), " is not a likley control and not in likelyMasks!"));
%         else
%             likelyMasks = obj.likelyMasks(num2str(u));
%         end
%     end
        
    for i = 1:length(obj.dims)
        sumPrior = 0.0;
        for j=1:length(obj.betas)-1
            sumPrior = sumPrior + x{j+2};
        end
        validProb = (sumPrior >= 0) .* (sumPrior <= 1);
        for j=1:length(obj.betas)-1
            validProb = validProb .* (x{j+2} >= 0) .* (x{j+2} <= 1);
        end
        if i == 1
            % NOTE: we need to "freeze" the dynamics when we have invalid
            % probabilities. 
            dx{i} = obj.v .* cos(u) .* validProb; %.* likelyMasks; % see how to incorporate quick sum across priors
        elseif i == 2
            dx{i} = obj.v .* sin(u) .* validProb; %.* likelyMasks;
        elseif i == 3
            betaPosteriors = obj.betaPosterior(x, u);
            % NOTE: These dynamics are for a STATIONARY beta model. 
            % TODO: check how to do beta update for multiple beta
            for j=1:length(obj.betas)-1
                dx{j+2} = obj.gamma * (betaPosteriors{j} - x{j+2}) .* ... 
                    validProb;
            end
%             if strcmp(obj.betaModel, 'static')
%                 % NOTE: These dynamics are for a STATIONARY beta model. 
%                 % TODO: check how to do beta update for multiple beta
%                 dx{i} = obj.gamma * (obj.betaPosterior(x, u, i - 2) - x{i}) .* ... 
%                     (x{i} >= 0) .* (x{i} <= 1); %.* likelyMasks;
%             else
%                 % NOTE: These dynamics are for a DYNAMIC beta model. 
%                 dx{i} = obj.gamma * ((obj.betaPosterior(x, u, i - 2) - x{i}) + ...
%                     (obj.alpha * x{i} + (1 - obj.alpha) * obj.DeltaB0 - x{i})).* ...
%                     (x{i} >= 0) .* (x{u} <= 1); %.* likelyMasks;
%             end  
        end
        
    end

    if returnVector
      dx = cell2mat(dx);
    end
end