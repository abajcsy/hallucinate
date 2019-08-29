classdef XYTTrajectory < handle
    properties
        contStates % Continuous (x, y, t) states.
        figh % Figure handle.
    end

    methods
        function obj = XYTTrajectory(contStates)
            obj.contStates = contStates;
            obj.figh = [];
        end

        function empty = isEmpty(obj)
            empty = (length(obj.contStates) == 0);
        end

        function valid = inBounds(obj, t)
            valid = ~obj.isEmpty() && ...
                    obj.contStates{1}(3) <= t && ...
                    t <= obj.contStates{end}(3);
        end

        function state = getState(obj, t)
            idx = 0;
            for i = 1:(length(obj.contStates) - 1)
                if t >= obj.contStates{i}(3) && t <= obj.contStates{i + 1}(3)
                    idx = i;
                    break;
                end
            end

            if idx == 0
                fprintf('(getState) Warning: time %f is out of bounds of trajectory!\n', ...
                        t);
                state = [0; 0; 0];
            else
                % Linearly interpolate between the states.
                s = (t - obj.contStates{idx}(3)) / ...
                    (obj.contStates{idx + 1}(3) - obj.contStates{idx}(3));
                state = obj.contStates{idx} + s * ...
                        (obj.contStates{idx + 1} - obj.contStates{idx});
            end
        end

        function h = draw(obj, showStates)
            h = {};
            for idx = 1:(length(obj.contStates) - 1)
                s = obj.contStates{idx};
                sNext = obj.contStates{idx + 1};
                h{end + 1} = plot([s(1), sNext(1)], ...
                                  [s(2), sNext(2)], ...
                                  'Color', [255, 148, 148]/255., 'Linewidth', ...
                                  2);
            end
        end
    end
end