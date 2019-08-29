classdef XYTState < handle
    properties
        x
        y
        t

        costToCome
        evalFunc
        parent
        closed
    end

    methods
        function obj = XYTState(x, y, t)
            obj.x = x;
            obj.y = y;
            obj.t = t;

            obj.costToCome = 0;
            obj.evalFunc = 0;
            obj.parent = 0;
            obj.closed = 0;
        end
    end
end
