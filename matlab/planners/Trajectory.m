classdef Trajectory < handle
    properties
        contStates % Continuous (x, y, theta, v, t) states.
        splines % Spline parameters (a, b, c, d).
    end

    methods
        function obj = Trajectory(contStates)
            obj.contStates = contStates;
            obj.splines = {};

            % Generate splines between each continuous state.
            for idx = 1:(length(contStates) - 1)
                [a, b, c, d] = unicycleThirdOrderSpline(contStates{idx}(1), ...
                                                        contStates{idx}(2), ...
                                                        contStates{idx}(3), ...
                                                        contStates{idx}(4), ...
                                                        contStates{idx + 1}(1), ...
                                                        contStates{idx + 1}(2), ...
                                                        contStates{idx + 1}(3), ...
                                                        contStates{idx + 1}(4), ...
                                                        contStates{idx + 1}(5) ...
                                                        - contStates{idx}(5));

                obj.splines{idx} = [a, b, c, d];
            end
        end

        function state = getState(obj, t)
            splineIdx = 0;
            tLow = 0;

            for idx = 1:(length(obj.contStates) - 1)
                if obj.contStates{idx}(5) <= t && t <= obj.contStates{idx + ...
                                        1}(5)
                    splineIdx = idx;
                    tLow = obj.contStates{idx}(5);
                    break;
                end
            end

            if splineIdx > 0
                p = obj.splines{splineIdx};
                s = t - tLow;

                x = p(1, 1) * s^3 + p(1, 2) * s^2 + p(1, 3) * s + p(1, 4);
                y = p(2, 1) * s^3 + p(2, 2) * s^2 + p(2, 3) * s + p(2, 4);
                theta = atan(y / x);

                vx = 3 * p(1, 1) * s^2 + 2 * p(1, 2) * s + p(1, 3);
                vy = 3 * p(2, 1) * s^2 + 2 * p(2, 2) * s + p(2, 3);
                v = norm([vx; vy]);

                state = [x, y, theta, v];
            else
                fprintf('Warning: time %f is out of bounds of trajectory!\n', ...
                        t);
                state = [0; 0; 0; 0];
            end
        end

        function control = getControl(obj, t)
            control = [0; 0];  % TODO Implement this
        end
    end
end