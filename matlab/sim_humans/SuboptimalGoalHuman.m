classdef SuboptimalGoalHuman < SimHuman
    %SUBOPTIMALGOALHUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xcurr
        uRange
        trueGoal
        waypts % cell arr
        currGoal
        currGoalIdx
    end
    
    methods
        function obj = SuboptimalGoalHuman(x0, v, uRange, trueGoal, waypts)
            %GAUSSIANGOALHUMAN Construct an instance of this class
            %   Detailed explanation goes here
            obj.xcurr = x0;
            obj.v = v;
            obj.uRange = uRange;
            obj.trueGoal = trueGoal;
            obj.waypts = waypts; 
            obj.currGoal = obj.waypts{1};
            obj.currGoalIdx = 1;
        end
        
        %% Simulates an action of the human according to the model.
        function [xnext, u] = simulate(obj, x, t, dt)
            eps = 0.05;
            if norm(x - obj.currGoal) < eps
                obj.currGoalIdx = obj.currGoalIdx + 1;
                if obj.currGoalIdx > length(obj.waypts)
                    obj.currGoalIdx = length(obj.waypts);
                end
                obj.currGoal = obj.waypts{obj.currGoalIdx};
            end
            
            u = atan2(obj.currGoal(2) - x(2), obj.currGoal(1) - x(1));
            
            % Figure out where the human would move in one timestep.
            xdot = obj.dynamics(obj.xcurr, u);
            xnext = obj.xcurr + dt*xdot;
            obj.xcurr = xnext;
        end
    end
end

