classdef LatticePlannerState < handle
  properties
    x
    y
    theta
    v
    t

    costToCome
    evalFunc
    parent
    closed
  end

  methods
    function obj = LatticePlannerState(x, y, theta, v, t)
      obj.x = x;
      obj.y = y;
      obj.theta = theta;
      obj.v = v;
      obj.t = t;

      obj.costToCome = 0;
      obj.evalFunc = 0;
      obj.parent = 0;
      obj.closed = 0;
    end
  end
end
