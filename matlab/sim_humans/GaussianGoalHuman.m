classdef GaussianGoalHuman < SimHuman
    %GAUSSIANGOALHUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xcurr
        sigma
        uRange
        trueGoal
    end
    
    methods
        function obj = GaussianGoalHuman(x0, v, sigma, uRange, trueGoal)
            %GAUSSIANGOALHUMAN Construct an instance of this class
            %   Detailed explanation goes here
            obj.xcurr = x0;
            obj.v = v;
            obj.sigma = sigma;
            obj.uRange = uRange;
            obj.trueGoal = trueGoal;
            
        end
        
        %% Simulates an action of the human according to the model.
        function [xnext, u] = simulate(obj, x, t, dt)
            
            uopt = atan2(obj.trueGoal{1}(2) - x(2), obj.trueGoal{1}(1) - x(1));
            
            if obj.sigma == 0.0
                u = uopt;
            else
                pd = makedist('Normal','mu', uopt ,'sigma', obj.sigma);

                % Truncate normal distribution to have finite support that is 
                % within the valid control bounds.
                clippedNorm = truncate(pd, obj.uRange(1), obj.uRange(2));

                % Generate random action from the normal distribution. 
                u = random(clippedNorm, 1, 1);
            end
            
            % Figure out where the human would move in one timestep.
            xdot = obj.dynamics(obj.xcurr, u);
            xnext = obj.xcurr + dt*xdot;
            obj.xcurr = xnext;
        end
    end
end

