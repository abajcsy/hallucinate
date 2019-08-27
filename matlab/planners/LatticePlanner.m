classdef LatticePlanner
    properties
        xDisc     % Discretization in x (m / unit).
        yDisc     % Discretization in y (m / unit).
        thetaDisc % Discretization in theta (rad / unit).
        vDisc     % Discretization in v (m / s / unit).
        tDisc     % Discretization in t (s / unit).

        staticObsMap % Occupancy grid of the static obstacles.
        dynObsMap % Time-dependent occupancy grid of the dynamic obstacles.

        states % Map of states.
        stateBounds % Bounds (i.e. box constraints) on the states.
        heurWeight % Heuristic weight / suboptimality bound for search.
    end

    methods
        function obj = LatticePlanner(xDisc, yDisc, thetaDisc, vDisc, tDisc, heurWeight)
            obj.xDisc = xDisc;
            obj.yDisc = yDisc;
            obj.thetaDisc = thetaDisc;
            obj.vDisc = vDisc;
            obj.tDisc = tDisc;

            obj.states = containers.Map();
            obj.stateBounds = containers.Map();
            obj.heurWeight = heurWeight;
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
            fprintf('Planning from start state (%f, %f, %f, %f, %f)\n', ...
                    startStateCont(1), ...
                    startStateCont(2), ...
                    startStateCont(3), ...
                    startStateCont(4), ...
                    startStateCont(5));
            % TODO (HACK) Need to clear the states.
            obj.states = containers.Map();

            startStateDisc = [obj.contToDisc(startStateCont(1), 1);
                              obj.contToDisc(startStateCont(2), 2);
                              wrapToPi(obj.contToDisc(startStateCont(3), 3));
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
            foundPath = 0;

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
                    path = obj.reconstructPath(state);

                    contStates = {};
                    for idx = 1:length(path)
                        contStates{idx} = [obj.discToCont(path{idx}.x, 1);
                                           obj.discToCont(path{idx}.y, 2);
                                           obj.discToCont(path{idx}.theta, 3);
                                           obj.discToCont(path{idx}.v, 4);
                                           obj.discToCont(path{idx}.t, 5)];
                        fprintf('State at t: %f is (%f, %f, %f, %f)\n', ...
                                contStates{idx}(5), contStates{idx}(1), ...
                                contStates{idx}(2), contStates{idx}(3), ...
                                contStates{idx}(4));
                    end
                    traj = Trajectory(contStates);
                    foundPath = 1;
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
                        succ.evalFunc = succ.costToCome + obj.heurWeight * ...
                            obj.heuristicCostToGoal(succ, goalXY);
                    end

                    if ~succ.closed
                        openList{length(openList) + 1} = succ;
                    end
                end
            end

            % If planner fails, return an empty trajectory.
            if ~foundPath
                fprintf('Planner failed to find a path after %d expansions!\n', ...
                        expansions);
                contStates = {};
                traj = Trajectory(contStates);
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

            % TODO These should be planner parameters.
            % localXRange = [1, 2];
            localXRange = [1, 3];
            localYRange = [-2, 2];
            aRange = [-1, 1];

            % TODO These should be parameters.
            f0 = 1;
            f1 = 1;

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
                                            obj.discToCont(a, 4)];

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
                        [a, b, c, d] = unicycleThirdOrderTimeSpline(startCoords(1), ...
                                                                    startCoords(2), ...
                                                                    startCoords(3), ...
                                                                    startCoords(4), ...
                                                                    f0, ...
                                                                    endCoords(1), ...
                                                                    endCoords(2), ...
                                                                    endCoords(3), ...
                                                                    endCoords(4), ...
                                                                    f1, ...
                                                                    obj ...
                                                                    .discToCont(deltaT, 5));

                        edges{idx} = [a, b, c, d];

                        % TODO Should be a planner parameter.
                        numSamples = 50;
                        if ~obj.isEdgeSafe(a, b, c, d, ...
                                           obj.discToCont(state.t, 5), ...
                                           obj.discToCont(state.t + deltaT, 5), ...
                                           numSamples)
                            continue;
                        end

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
            valid = valid && obj.inBounds(t, 't');
        end

        function valid = inBounds(obj, value, stateName)
            valid = 1;

            if obj.stateBounds.isKey(stateName)
                bounds = obj.stateBounds(stateName);
                if value < bounds(1) || value > bounds(2)
                    valid = 0;
                end
            else
                fprintf('Warning: no bounds found for state name %s!\n', ...
                        stateName);
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
            cost = 0;

            % Add a cost on distance.
            cost = cost + norm([obj.discToCont(succ.x, 1); obj.discToCont(succ.y, 2)] - ...
                               [obj.discToCont(state.x, 1); obj.discToCont(state.y, ...
                                                              2)]);

            % TODO Add additional costs as needed
        end

        function safe = isEdgeSafe(obj, a, b, c, d, t0, t1, numSamples)
            safe = 1;

            pfunc = @(t) a(3) .* t.^3 + b(3) .* t.^2 + c(3) .* t + d(3);

            xfunc = @(t) a(1) .* pfunc(t).^3 + b(1) .* pfunc(t).^2 + c(1) .* ...
                    pfunc(t) + d(1);
            yfunc = @(t) a(2) .* pfunc(t).^3 + b(2) .* pfunc(t).^2 + c(2) .* ...
                    pfunc(t) + d(2);

            T = t1 - t0;
            Tsample = T / numSamples;
            for t = 1:(numSamples - 1)
                xsample = xfunc(t * Tsample);
                ysample = yfunc(t * Tsample);

                % Check against the static obstacle map.
                if obj.staticObsMap.getData(xsample, ysample, 0) > 0
                    safe = 0;
                    break;
                end

                % Check against the dynamic obstacle map.
                if obj.dynObsMap.getData(xsample, ysample, t0 + Tsample * t) ...
                        > 0
                    safe = 0;
                    break;
                end
            end
        end

        function isVelocityBounded(obj, a, b, c, d, t0, t1, numSamples)
        % TODO Finish implementing / testing this function.
            safe = 1;
            vxfunc = @(t) 3 .* a(1) .* t.^2 + 2 .* b(1) .* t + ...
                     c(1);
            vyfunc = @(t) 3 .* a(2) .* t.^2 + 2 .* b(2) .* t + ...
                     c(2);

            T = t1 - t0;
            Tsample = T / numSamples;
            for t = 1:(numSamples - 1)
                vxsample = vxfunc(t * Tsample);
                vysample = vyfunc(t * Tsample);

                % Check the sampled velocity to ensure that it is in bounds.
                vsample = sqrt(vxsample^2 + vysample^2);
                if ~obj.inBounds(vsample, 'v')
                    safe = 0;
                    break;
                end
            end
        end
    end
end