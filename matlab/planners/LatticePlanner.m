classdef LatticePlanner < PlannerBase
  properties
    xDisc
    yDisc
    thetaDisc
    vDisc
    tDisc
  end

  enumeration
    STATE_X
    STATE_Y
    STATE_THETA
    STATE_V
    STATE_T
  end

  methods
    function obj = LatticePlanner(env)
      obj@PlannerBase(env);
    end

    function traj = plan()
      % TODO
      traj = [];
    end

    function contVal = discToCont(obj, discVal, var)
      switch var
        case LatticePlanner.STATE_X
          contVal = obj.xDisc * discVal;
        case LatticePlanner.STATE_Y
          contVal = obj.yDisc * discVal;
        case LatticePlanner.STATE_THETA
          contVal = obj.thetaDisc * discVal;
        case LatticePlanner.STATE_V
          contVal = obj.vDisc * discVal;
        case LatticePlanner.STATE_T
          contVal = obj.tDisc * discVal;
        otherwise
          fprintf("Error: no such state variable %d", var);
      end
    end

    % TODO discVal = contToDisc(obj, contVal, var) ...
    % TODO cost = getCost(a, b, c, d)
  end
end
