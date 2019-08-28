classdef FixedTrajHuman < SimHuman
    %FIXEDTRAJHUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0  % initial state
        ctrls    % controls
        times   % time interval in which to apply each control
    end
    
    methods
        function obj = FixedTrajHuman(x0, v, ctrls, times)
            %FIXEDTRAJHUMAN Construct an instance of this class
            %   Detailed explanation goes here
            obj.x0 = x0;
            obj.v = v;
            obj.ctrls = ctrls;
            obj.times = times;
        end
        
        %% Simulate the human action. 
        function [xnext, u] = simulate(obj, x, t, dt)
            u = nan;
            for i=1:length(obj.times)
                min = obj.times{i}(1);
                max = obj.times{i}(2);
                if t >= min && t < max
                    u = obj.ctrls{i};
                    break;
                end
            end
            
            % Time is out of bounds. 
            if isnan(u)
                minL = obj.times{1}(1);
                maxU = obj.times{end}(2);
                if t < minL
                    u = obj.ctrls{1};
                elseif t >= maxU 
                    u = obj.ctrls{end};
                else
                    warning('Time %f does not make sense!\n', t);
                    u = 0;
                end
            end
            
            % Figure out where the human would move in one timestep.
            xdot = obj.dynamics(x, u);
            xnext = x + dt*xdot;
        end

    end
end

