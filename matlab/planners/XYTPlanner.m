classdef XYTPlanner
    properties
        xDisc     % Discretization in x (m / unit).
        yDisc     % Discretization in y (m / unit).
        tDisc     % Discretization in t (s / unit).
        deltaT    % Number of discrete time steps between a state and it's successor.

        staticObsMap % Occupancy grid of the static obstacles.
        dynObsMap % Time-dependent occupancy grid of the dynamic obstacles.

        states % Map of states.
        stateBounds % Bounds (i.e. box constraints) on the states.
        heurWeight % Heuristic weight / suboptimality bound for search.
    end

    methods
        function obj = XYTPlanner(xDisc, yDisc, tDisc, heurWeight, deltaTCont)
            obj.xDisc = xDisc;
            obj.yDisc = yDisc;
            obj.tDisc = tDisc;

            obj.states = containers.Map();
            obj.stateBounds = containers.Map();

            obj.heurWeight = heurWeight;
            obj.deltaT = obj.contToDisc(deltaTCont, 3);
            fprintf('%f s (%d units)\n', deltaTCont, obj.deltaT);

            obj.staticObsMap = 0;
            obj.dynObsMap = 0;
        end

        function state = getState(obj, xDisc, yDisc, tDisc)
            key = mat2str([xDisc, yDisc, tDisc]);
            if obj.states.isKey(key)
                state = obj.states(key);
            else
                state = XYTState(xDisc, yDisc, tDisc);
                state.costToCome = 1e9;
                obj.states(key) = state;
            end
        end

        function traj = plan(obj, startStateCont, goalXY, goalTol)
            fprintf('Planning from start state (%f, %f, %f)\n', ...
                    startStateCont(1), ...
                    startStateCont(2), ...
                    startStateCont(3));
            % TODO (HACK) Need to clear the states.
            obj.states = containers.Map();

            startStateDisc = [obj.contToDisc(startStateCont(1), 1);
                              obj.contToDisc(startStateCont(2), 2);
                              obj.contToDisc(startStateCont(3), 3)];

            % TODO Priority queue???
            startState = obj.getState(startStateDisc(1), ...
                                      startStateDisc(2), ...
                                      startStateDisc(3));
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
                                           obj.discToCont(path{idx}.t, 3)];
                        fprintf('State at t: %f is (%f, %f)\n', ...
                                contStates{idx}(3), contStates{idx}(1), ...
                                contStates{idx}(2));
                    end
                    traj = XYTTrajectory(contStates);
                    foundPath = 1;
                    break;
                end

                % Otherwise, expand the state.
                [succs, costs] = obj.expand(state);
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
                traj = XYTTrajectory(contStates);
            end
        end

        function heur = heuristicCostToGoal(obj, state, goalXY)
            xCont = obj.discToCont(state.x, 1);
            yCont = obj.discToCont(state.y, 2);

            heur = norm(goalXY - [xCont; yCont]);
        end

        function [succs, costs] = expand(obj, state)
            succs = {};
            costs = {};

            idx = 1;

            for deltaX = -1:1:1
                for deltaY = -1:1:1
                    if deltaX == 0 && deltaY == 0
                        continue;
                    end

                    startCoords = [obj.discToCont(state.x, 1);
                                   obj.discToCont(state.y, 2);
                                   obj.discToCont(state.t, 3)];

                    endCoords = [obj.discToCont(state.x + deltaX, 1);
                                 obj.discToCont(state.y + deltaY, 2);
                                 obj.discToCont(state.t + obj.deltaT, 3)];

                    numSamples = 20;

                    % If the coordinates of the successor state are not
                    % valid, do not generate the edge.
                    if ~obj.isStateValid(endCoords)
                        continue;
                    end

                    if ~obj.isEdgeValid(startCoords, endCoords, numSamples)
                        continue;
                    end

                    % Add the successor state.
                    succ = obj.getState(obj.contToDisc(endCoords(1), 1), ...
                                        obj.contToDisc(endCoords(2), 2), ...
                                        obj.contToDisc(endCoords(3), 3));

                    costs{idx} = obj.getCost(state, succ);
                    succs{idx} = succ;

                    idx = idx + 1;
                end
            end
        end

        function valid = isStateValid(obj, stateCoords)
            valid = 1;

            % Check state bounds provided by the environment.
            valid = valid && obj.inBounds(stateCoords(1), 'x');
            valid = valid && obj.inBounds(stateCoords(2), 'y');
            valid = valid && obj.inBounds(stateCoords(3), 't');
        end

        function valid = isEdgeValid(obj, stateCoords, succCoords, ...
                                     numSamples)
            valid = 1;

            % Check against collisions with dynamic / static obstacles.
            totalDist = norm(succCoords(1:2) - stateCoords(1:2));
            sampleDist = totalDist / numSamples;
            for k = 0:(numSamples - 1)
                sampleCoords = stateCoords + (sampleDist * k) * (succCoords - ...
                                                                 stateCoords);

                % Check against the static obstacle map.
                if obj.staticObsMap ~= 0 && obj.staticObsMap.getData(sampleCoords(1), ...
                                                                     sampleCoords(2), 0) > 0
                    valid = 0;
                    break;
                end

                % Check against the dynamic obstacle map.
                if obj.dynObsMap ~= 0 && obj.dynObsMap.getData(sampleCoords(1), ...
                                                               sampleCoords(2), ...
                                                               sampleCoords(3)) ...
                        > 0
                    valid = 0;
                    break;
                end
            end
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
                discVal = contVal / obj.tDisc;
              otherwise
                fprintf("Error: no such state variable %d", var);
            end
        end

        function cost = getCost(obj, state, succ)
            cost = 0;

            % Add a cost on distance.
            cost = cost + norm([obj.discToCont(succ.x, 1); obj.discToCont(succ.y, 2)] - ...
                               [obj.discToCont(state.x, 1); obj.discToCont(state.y, ...
                                                              2)]);
        end
    end
end