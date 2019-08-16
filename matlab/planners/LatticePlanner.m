classdef LatticePlanner < PlannerBase
  properties
    xDisc
    yDisc
    thetaDisc
    vDisc
    tDisc
  end

  methods
    function obj = LatticePlanner(env, xDiscArg, yDiscArg, thetaDiscArg, vDiscArg, tDiscArg)
      obj@PlannerBase(env);

      obj.xDisc = xDiscArg;
      obj.yDisc = yDiscArg;
      obj.thetaDisc = thetaDiscArg;
      obj.vDisc = vDiscArg;
      obj.tDisc = tDiscArg;
    end

    function traj = plan(obj, startStateCont, goalXY, goalTol)
      startStateDisc = [obj.contToDisc(startStateCont(1), 1);
                        obj.contToDisc(startStateCont(2), 2);
                        obj.contToDisc(startStateCont(3), 3);
                        obj.contToDisc(startStateCont(4), 4);
                        obj.contToDisc(startStateCont(5), 5)];

      % TODO Priority queue???
      startState = LatticePlannerState(startStateDisc(1), ...
                                       startStateDisc(2), ...
                                       startStateDisc(3), ...
                                       startStateDisc(4), ...
                                       startStateDisc(5));
      startState.evalFunc = obj.heuristicCostToGoal(startState, goalXY);
      openList = {startState};

      expansions = 0;

      while length(openList) > 0
        % Search for the state with the smallest evaluation function.
        bestIdx = -1;
        bestEvalFunc = 1e9;
        for i = 1:length(openList)
          if openList{i}.evalFunc < bestEvalFunc
            bestIdx = i;
            bestEvalFunc = openList{i}.evalFunc;
          end
        end

        state = openList{bestIdx};
        openList(:, bestIdx) = [];

        if mod(expansions, 100) == 0
          fprintf("Expanded %d states so far\n", expansions);
        end

        % Check if the goal has been reached.
        contX = obj.discToCont(state.x, 1);
        contY = obj.discToCont(state.y, 2);
        if norm(goalXY - [contX; contY]) < goalTol
          fprintf("Found goal after %d expansions\n", expansions);
          traj = obj.reconstructPath(state)
          break;
        end

        % Otherwise, expand the state.
        [succs, edges, costs] = obj.expand(state);
        expansions = expansions + 1;
        for i = 1:length(succs)
          succ = succs{i};
          succ.costToCome = state.costToCome + costs{i};
          succ.parent = state;
          succ.evalFunc = succ.costToCome + obj.heuristicCostToGoal(succ, goalXY);
          openList{length(openList) + 1} = succ;
        end
      end
    end

    function heur = heuristicCostToGoal(obj, state, goalXY)
      xCont = obj.discToCont(state.x, 1);
      yCont = obj.discToCont(state.y, 2);

      heur = norm(goalXY - [xCont; yCont]);
    end

    function [succs, edges, costs] = expand(obj, state)
      succs = {};
      edges = {};
      costs = {};

      idx = 1;
      deltaT = 1;

      for deltaXDisc = 0:1:1
        for deltaYDisc = -1:1:1
          for deltaVDisc = -1:1:1
            for deltaThetaDisc = -1:1:1
              if deltaXDisc ~= 0 || deltaYDisc ~= 0
                succ = LatticePlannerState(state.x + deltaXDisc, ...
                                           state.y + deltaYDisc, ...
                                           state.theta + deltaThetaDisc, ...
                                           state.v + deltaVDisc, ...
                                           state.t + deltaT);
                if ~obj.isValidState(succ)
                  continue;
                end

                [a, b, c, d] = unicycleThirdOrderSpline(obj.discToCont(state.x, 1), ...
                                                        obj.discToCont(state.y, 2), ...
                                                        obj.discToCont(state.theta, 3), ...
                                                        obj.discToCont(state.v, 4), ...
                                                        obj.discToCont(succ.x, 1), ...
                                                        obj.discToCont(succ.y, 2), ...
                                                        obj.discToCont(succ.theta, 3), ...
                                                        obj.discToCont(succ.v, 4), ...
                                                        obj.discToCont(deltaT, 5));

                edges{idx} = [a, b, c, d];
                % costs{idx} = 1; % TODO
                % costs{idx} = abs(deltaXDisc) + abs(deltaYDisc) + abs(deltaThetaDisc) + abs(deltaVDisc);
                % costs{idx} = sqrt((obj.discToCont(succ.x, 1) - obj.discToCont(state.x, 1))^2 + ...
                %                   (obj.discToCont(succ.y, 2) - obj.discToCont(state.y, 2))^2);
                costs{idx} = norm([obj.discToCont(succ.x, 1); obj.discToCont(succ.y, 2)] - ...
                                  [obj.discToCont(state.x, 1); obj.discToCont(state.y, 2)]);
                succs{idx} = succ;

                idx = idx + 1;
              end
            end
          end
        end
      end
    end

    function valid = isValidState(obj, state)
      valid = 1;

      if state.v < 0
        valid = 0;
      end

      % TODO Here also check environment boundary.
    end

    function path = reconstructPath(obj, goalState)
      fprintf("Reconstructing the path...\n");
      path = {};
      idx = 1;

      p = goalState;
      while p ~= 0
        path{idx} = p;
        p = p.parent;
        idx = idx + 1;
      end

      path = flip(path);
    end

    function contVal = discToCont(obj, discVal, var)
      switch var
        case 1
          contVal = obj.xDisc * discVal;
        case 2
          contVal = obj.yDisc * discVal;
        case 3
          contVal = obj.thetaDisc * discVal;
        case 4
          contVal = obj.vDisc * discVal;
        case 5
          contVal = obj.tDisc * discVal;
        otherwise
          fprintf("Error: no such state variable %d", var);
      end
    end

    function discVal = contToDisc(obj, contVal, var)
      switch var
        case 1
          discVal = contVal / obj.xDisc;
        case 2
          discVal = contVal / obj.yDisc;
        case 3
          discVal = contVal / obj.thetaDisc;
        case 4
          discVal = contVal / obj.vDisc;
        case 5
          discVal = contVal / obj.tDisc;
        otherwise
          fprintf("Error: no such state variable %d", var);
      end
    end

    function cost = getCost(a, b, c, d)
      cost = 1;
    end
  end
end
