classdef LatticePlanner
    properties
        xDisc
        yDisc
        thetaDisc
        vDisc
        tDisc

        states
        stateBounds
    end

    methods
        function obj = LatticePlanner(xDisc, yDisc, thetaDisc, vDisc, tDisc)
            obj.xDisc = xDisc;
            obj.yDisc = yDisc;
            obj.thetaDisc = thetaDisc;
            obj.vDisc = vDisc;
            obj.tDisc = tDisc;

            obj.states = containers.Map();
            obj.stateBounds = containers.Map();
        end

        function state = getState(obj, xDisc, yDisc, thetaDisc, vDisc, tDisc)
            key = mat2str([xDisc, yDisc, thetaDisc, vDisc, tDisc]);
            if obj.states.isKey(key)
                state = obj.states(key);
            else
                state = LatticePlannerState(xDisc, yDisc, thetaDisc, vDisc, tDisc);
                state.costToCome = 1e9;
                obj.states(key) = state;
            end
        end

        function traj = plan(obj, startStateCont, goalXY, goalTol)
            startStateDisc = [obj.contToDisc(startStateCont(1), 1);
                              obj.contToDisc(startStateCont(2), 2);
                              obj.contToDisc(startStateCont(3), 3);
                              obj.contToDisc(startStateCont(4), 4);
                              obj.contToDisc(startStateCont(5), 5)];

            % TODO Priority queue???
            startState = obj.getState(startStateDisc(1), ...
                                      startStateDisc(2), ...
                                      startStateDisc(3), ...
                                      startStateDisc(4), ...
                                      startStateDisc(5));
            startState.costToCome = 0;
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
                if state.closed
                    continue;
                end

                if mod(expansions, 100) == 0
                    fprintf("Expanded %d states so far\n", expansions);
                end

                % Check if the goal has been reached.
                contX = obj.discToCont(state.x, 1);
                contY = obj.discToCont(state.y, 2);
                % norm(goalXY - [contX; contY])
                if norm(goalXY - [contX; contY]) < goalTol
                    fprintf("Found goal after %d expansions\n", expansions);
                    traj = obj.reconstructPath(state);
                    break;
                end

                % Otherwise, expand the state.
                [succs, edges, costs] = obj.expand(state);
                % TODO visualize the expansion
                state.closed = 1;
                expansions = expansions + 1;
                for i = 1:length(succs)
                    succ = succs{i};
                    if state.costToCome + costs{i} < succ.costToCome
                        succ.costToCome = state.costToCome + costs{i};
                        succ.parent = state;
                        succ.evalFunc = succ.costToCome + obj.heuristicCostToGoal(succ, goalXY);
                    end

                    if ~succ.closed
                        openList{length(openList) + 1} = succ;
                    end
                end
            end
        end

        function heur = heuristicCostToGoal(obj, state, goalXY)
            xCont = obj.discToCont(state.x, 1);
            yCont = obj.discToCont(state.y, 2);

            heur = norm(goalXY - [xCont; yCont]);
            % heur = 0;
        end

        function [succs, edges, costs] = expand(obj, state)
            succs = {};
            edges = {};
            costs = {};

            idx = 1;
            deltaT = 1;

            % localXRange = [1, 2];
            localXRange = [1, 3];
            localYRange = [-2, 2];
            aRange = [-1, 1];

            for localX = localXRange(1):localXRange(2)
                for localY = localYRange(1):localYRange(2)
                    for a = aRange(1):aRange(2)
                        % Get the state and successor coordinates in the
                        % global frame.
                        startCoords = [obj.discToCont(state.x, 1);
                                       obj.discToCont(state.y, 2);
                                       obj.discToCont(state.theta, 3);
                                       obj.discToCont(state.v, 4)];

                        R = [cos(startCoords(3)), -sin(startCoords(3));
                             sin(startCoords(3)), cos(startCoords(3))];

                        localXCont = obj.discToCont(localX, 1);
                        localYCont = obj.discToCont(localY, 2);

                        endCoords = startCoords + [R * [localXCont; localYCont]; ...
                                            atan(localYCont / localXCont); ...
                                            a];

                        % If the coordinates of the successor state are not
                        % valid, do not generate the edge.
                        if ~obj.isValid(endCoords(1), ...
                                        endCoords(2), ...
                                        endCoords(3), ...
                                        endCoords(4), ...
                                        obj.discToCont(state.t + deltaT, 5))
                            continue;
                        end

                        % Otherwise, generate motion primitive.
                        [a, b, c, d] = unicycleThirdOrderSpline(startCoords(1), ...
                                                                startCoords(2), ...
                                                                startCoords(3), ...
                                                                startCoords(4), ...
                                                                endCoords(1), ...
                                                                endCoords(2), ...
                                                                endCoords(3), ...
                                                                endCoords(4), ...
                                                                obj ...
                                                                .discToCont(deltaT, 5));

                        edges{idx} = [a, b, c, d];

                        % Add the successor state.
                        succ = obj.getState(obj.contToDisc(endCoords(1), 1), ...
                                            obj.contToDisc(endCoords(2), 2), ...
                                            obj.contToDisc(endCoords(3), 3), ...
                                            obj.contToDisc(endCoords(4), 4), ...
                                            state.t + deltaT);

                        costs{idx} = obj.getCost(state, succ, a, b, ...
                                                        c, d);
                        succs{idx} = succ;

                        idx = idx + 1;
                    end
                end
            end
        end

        function valid = isValid(obj, x, y, theta, v, t)
            valid = 1;

            % Check state bounds provided by the environment.
            valid = valid && obj.inBounds(x, 'x');
            valid = valid && obj.inBounds(y, 'y');
            valid = valid && obj.inBounds(theta, 'theta');
            valid = valid && obj.inBounds(v, 'v');
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

        function cost = getCost(obj, state, succ, a, b, c, d)
            cost = norm([obj.discToCont(succ.x, 1); obj.discToCont(succ.y, 2)] - ...
                        [obj.discToCont(state.x, 1); obj.discToCont(state.y, ...
                                                              2)]);

            numSamples = 50;
            T = obj.discToCont(succ.t - state.t, 5);
            % cost = cost + 10 * obj.getSteeringCost(a, b, c, d, T, numSamples);
        end

        function cost = getSteeringCost(obj, a, b, c, d, T, numSamples)
          xfunc = @(t) a(1) .* t.^3 + b(1) .* t.^2 + ...
                  c(1) .* t + d(1);
          yfunc = @(t) a(2) .* t.^3 + b(2) .* t.^2 + ...
                  c(2) .* t + d(2);
          thetafunc = @(t) atan(yfunc(t) / xfunc(t));

          Tsample = T / numSamples;
          cost = 0;
          for t = 1:(numSamples - 1)
              cost = cost + abs(thetafunc(t + 1) - thetafunc(t)) / Tsample;
          end
        end
    end
end