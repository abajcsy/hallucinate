classdef Environment
    properties
        stateBounds % Bounds on the states.
    end

    methods
        function obj = Environment()
            obj.stateBounds = containers.Map();
        end

        function valid = inBounds(obj, value, stateName)
            valid = 1;

            if obj.stateBounds.isKey(stateName)
                bounds = obj.stateBounds(stateName);
                if value < bounds(1) || value > bounds(2)
                    valid = 0;
                end
            end
        end
    end
end
