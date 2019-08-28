classdef SimHuman < handle
    %SIMHUMAN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        v
    end
    
    methods

        % No constructor in DynSys class. Use constructors in the subclasses
        
        %% Simulates the next state and action. 
        function [xnext, u] = simulate(obj, x, t, dt)
            error("Not implemented in base class!");
        end
        
        %% Defines the dynamics of the human. 
        function xdot = dynamics(obj, x, u)
            % Simple dynamics of human.
            xdot = [obj.v*cos(u); obj.v*sin(u)];
        end
    end
end

