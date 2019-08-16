classdef (Abstract) PlannerBase < handle
  properties
    env
  end

  methods
    function obj = PlannerBase(env)
      obj.env = env;
    end
  end

  methods (Abstract)
    traj = plan()
  end
end
