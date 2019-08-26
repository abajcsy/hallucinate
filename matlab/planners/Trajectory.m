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

                vx = 3 * p(1, 1) * s^2 + 2 * p(1, 2) * s + p(1, 3);
                vy = 3 * p(2, 1) * s^2 + 2 * p(2, 2) * s + p(2, 3);
                v = norm([vx; vy]);

                theta = atan2(vy, vx);

                state = [x; y; theta; v];
            else
                fprintf('Warning: time %f is out of bounds of trajectory!\n', ...
                        t);
                state = [0; 0; 0; 0];
            end
        end

        function control = getControl(obj, t)
            splineIdx = 0;
            tLow = 0;

            omega = 0;
            a = 0;

            for idx = 1:(length(obj.contStates) - 1)
                if obj.contStates{idx}(5) <= t && t <= obj.contStates{idx + ...
                                        1}(5)
                    splineIdx = idx;
                    tLow = obj.contStates{idx}(5);
                    break;
                end
            end

            fprintf('At t = %f, tLow = %f, on spline index = %d\n', ...
                    t, tLow, splineIdx);

            if splineIdx > 0
                p = obj.splines{splineIdx};
                s = t - tLow;

                a = [p(1, 1); p(2, 1)];
                b = [p(1, 2); p(2, 2)];
                c = [p(1, 3); p(2, 3)];
                d = [p(1, 4); p(2, 4)];

                vx = 3 * a(1) * s^2 + 2 * b(1) * s + c(1);
                vy = 3 * a(2) * s^2 + 2 * b(2) * s + c(2);

                v = norm([vx; vy]);
                theta = atan2(vy, vx);

                omega = ((6 * a(2) * s + 2 * b(2)) * cos(theta) - ...
                         (6 * a(1) * s + 2 * b(1)) * sin(theta)) / v;
                a = (v * omega * sin(theta) + 6 * a(1) * s + 2 * b(1)) ...
                    / cos(theta);
            else
                fprintf('Warning: time %f is out of bounds of trajectory!\n', ...
                        t);
            end

            % Return the control input.
            control = [omega; a];
        end

        function draw(obj, showStates)
            for idx = 1:length(obj.splines)
                s = obj.contStates{idx};
                sNext = obj.contStates{idx + 1};
                tDisc = sNext(5) - s(5);

                if nargin > 1 && showStates
                    obj.drawTriangle([s(1); s(2)], s(3), 0.1);
                end

                p = obj.splines{idx};

                a = [p(1, 1); p(2, 1)];
                b = [p(1, 2); p(2, 2)];
                c = [p(1, 3); p(2, 3)];
                d = [p(1, 4); p(2, 4)];

                xfunc = @(t) a(1) .* t.^3 + b(1) .* t.^2 + c(1) .* t + d(1);
                yfunc = @(t) a(2) .* t.^3 + b(2) .* t.^2 + c(2) .* t + d(2);

                fplot(xfunc, yfunc, [0 tDisc]);
            end
        end

        function tri = drawTriangle(obj, origin, rot, sideLength)
            aLocal = [2*sideLength; 0; 1];
            bLocal = [0; -sideLength / sqrt(2); 1];
            cLocal = [0; sideLength / sqrt(2); 1];

            T = [cos(rot), -sin(rot), origin(1);
                 sin(rot), cos(rot), origin(2);
                 0, 0, 1];
            a = T * aLocal;
            b = T * bLocal;
            c = T * cLocal;

            tri = fill([a(1), b(1), c(1)], ...
                       [a(2), b(2), c(2)], ...
                       'b');
            set(tri, 'facealpha', 0.5);
        end
    end
end