classdef Unicycle4DRobot < DynSys
    %DUBINSCARROBOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Angular control bounds
        wMax

        % Acceleration control bounds
        aRange

        % Speed control bounds
        vRange
        
        dims
    end
    
    methods
        function obj = Unicycle4DRobot(x, wMax, aRange, vRange)
            % obj = Unicycle4DRobot(x, wMax, aRange, dMax)
            %     Dynamics of the Unicycle4DRobot
            %         \dot{x}_1 = x_4 * cos(x_3) 
            %         \dot{x}_2 = x_4 * sin(x_3) 
            %         \dot{x}_3 = u_1 = w
            %         \dot{x}_4 = u_2 = a
            %           wMin <= w <= wMax
            %           aMin <= a <= aMax
            obj.wMax = wMax;
            obj.aRange = aRange;
            obj.vRange = vRange;
            
            obj.x = x;
            
            obj.nx = length(x);
            obj.nu = 2;
            obj.dims = 1:4;
        end
    end
end

