classdef Trajectory < handle
    properties
        contStates % Continuous (x, y, theta, v, t) states.
        splines % Spline parameters (a, b, c, d).
        figh % figure handle
    end

    methods
        function obj = Trajectory(contStates)
            obj.contStates = contStates;
            obj.splines = {};
            obj.figh = [];

            % Generate splines between each continuous state.
            for idx = 1:(length(contStates) - 1)
                [f0, f1] = obj.getHeuristicFValues(contStates{idx}(4), ...
                                                   contStates{idx}(1), ...
                                                   contStates{idx}(2), ...
                                                   contStates{idx + 1}(1), ...
                                                   contStates{idx + 1}(2));

                [a, b, c, d] = unicycleThirdOrderTimeSpline(contStates{idx}(1), ...
                                                            contStates{idx}(2), ...
                                                            contStates{idx}(3), ...
                                                            contStates{idx}(4), ...
                                                            f0, ...
                                                            contStates{idx + 1}(1), ...
                                                            contStates{idx + 1}(2), ...
                                                            contStates{idx + 1}(3), ...
                                                            contStates{idx + 1}(4), ...
                                                            f1, ...
                                                            contStates{idx + 1}(5) ...
                                                            - contStates{idx}(5));

                obj.splines{idx} = [a, b, c, d];
            end
        end

        function empty = isEmpty(obj)
            empty = isempty(obj.splines);
        end

        function [f0, f1] = getHeuristicFValues(obj, v0, x0, y0, x1, y1)
            f0 = v0 + norm([x1 - x0; y1 - y0]);
            f1 = v0 + norm([x1 - x0; y1 - y0]);
        end
        
        function splineIdx = splineIndex(obj, t)
            splineIdx = 0;

            for idx = 1:(length(obj.contStates) - 1)
                if obj.contStates{idx}(5) <= t && t <= obj.contStates{idx + ...
                                        1}(5)
                    splineIdx = idx;
                    break;
                end
            end
        end

        function [a, b, c, d, tLow] = getSpline(obj, t)
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

            if splineIdx <= 0
                fprintf('Warning: time %f is out of bounds of trajectory!\n', ...
                        t);
                a = zeros(3, 1);
                b = zeros(3, 1);
                c = zeros(3, 1);
                d = zeros(3, 1);
            else
                p = obj.splines{splineIdx};

                a = [p(1, 1); p(2, 1); p(3, 1)];
                b = [p(1, 2); p(2, 2); p(3, 2)];
                c = [p(1, 3); p(2, 3); p(3, 3)];
                d = [p(1, 4); p(2, 4); p(3, 4)];
            end
        end

        function state = getState(obj, t)
            [a, b, c, d, tLow] = obj.getSpline(t);

            % Time always starts from zero.
            s = t - tLow;

            p = a(3) * s^3 + b(3) * s^2 + c(3) * s + d(3);
            x = a(1) * p^3 + b(1) * p^2 + c(1) * p + d(1);
            y = a(2) * p^3 + b(2) * p^2 + c(2) * p + d(2);

            vp = 3 * a(3) * s^2 + 2 * b(3) * s + c(3);
            vx = (3 * a(1) * p^2 + 2 * b(1) * p + c(1)) * vp;
            vy = (3 * a(2) * p^2 + 2 * b(2) * p + c(2)) * vp;
            v = norm([vx; vy]);

            theta = atan2(vy, vx);

            state = [x; y; theta; v];
        end

        function control = getControl(obj, t)
            [a, b, c, d, tLow] = obj.getSpline(t);

            % Time always starts from zero.
            s = t - tLow;

            omega = 0;
            acc = 0;

            p = a(3) * s^3 + b(3) * s^2 + c(3) * s + d(3);
            vp = 3 * a(3) * s^2 + 2 * b(3) * s + c(3);
            ap = 6 * a(3) * s + 2 * b(3);

            vxp = (3 * a(1) * p^2 + 2 * b(1) * p + c(1));
            vyp = (3 * a(2) * p^2 + 2 * b(2) * p + c(2));
            axp = (6 * a(1) * p + 2 * b(1));
            ayp = (6 * a(2) * p + 2 * b(2));
            vx = vxp * vp;
            vy = vyp * vp;
            v = norm([vx; vy]);

            ax = (6 * a(1) * p + 2 * b(1)) * vp^2 + (3 * a(1) * p^2 + 2 * b(1) ...
                                                     * p + c(1)) * ap;
            ay = (6 * a(2) * p + 2 * b(2)) * vp^2 + (3 * a(2) * p^2 + 2 * b(2) ...
                                                     * p + c(2)) * ap;

            theta = atan2(vyp, vxp);

            % acc = ap * norm([vxp; vyp]) + (vp^2 / norm([vxp; vyp])) * (axp + ...
            %                                                   ayp);
            % omega = (vp / (vxp^2 + vyp^2)) * (ayp * vxp - axp * vyp);

            omega = (ay * cos(theta) - ax * sin(theta)) / v;
            acc = (ax + omega * v * sin(theta)) / cos(theta);

            % Return the control input.
            control = [omega; acc];
        end

        function control = getProportionalControl(obj, t, x, kPropAngAcc, ...
                                                       kPropLinAcc)
            xRef = obj.getState(t);
            posErr = [xRef(1) - x(1);
                      xRef(2) - x(2)];

            % Acceleration is proportional to position error along current
            % direction of travel.
            dir = [cos(x(3)); sin(x(3))];
            a = kPropLinAcc * (posErr' * dir);

            % Angular velocity is proportional to the difference between the
            % current heading and the reference heading.
            headingErr = xRef(3) - x(3);
            omega = kPropAngAcc * headingErr;

            control = [omega; a];
        end

        function h = draw(obj, showStates)
            h = {};
            for idx = 1:length(obj.splines)
                s = obj.contStates{idx};
                sNext = obj.contStates{idx + 1};
                tDisc = sNext(5) - s(5);

                if nargin > 1 && showStates
                    obj.drawTriangle([s(1); s(2)], s(3), 0.1);
                end

                p = obj.splines{idx};

                a = [p(1, 1); p(2, 1); p(3, 1)];
                b = [p(1, 2); p(2, 2); p(3, 2)];
                c = [p(1, 3); p(2, 3); p(3, 3)];
                d = [p(1, 4); p(2, 4); p(3, 4)];

                pfunc = @(t) a(3) .* t.^3 + b(3) .* t.^2 + c(3) .* t + d(3);
                xfunc = @(t) a(1) .* pfunc(t).^3 + b(1) .* pfunc(t).^2 + c(1) .* pfunc(t) + d(1);
                yfunc = @(t) a(2) .* pfunc(t).^3 + b(2) .* pfunc(t).^2 + c(2) .* pfunc(t) + d(2);

                h{end+1} = fplot(xfunc, yfunc, [0 tDisc], ...
                    'Color', [255, 148, 148]/255., 'Linewidth', 2);
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