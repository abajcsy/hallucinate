function [a, b, c, d] = unicycleThirdOrderTimeSpline(x0, y0, th0, v0, f0, x1, ...
                                                     y1, th1, v1, f1, T)
%     fprintf("Computing spline from (%f, %f, %f, %f) to (%f, %f, %f, %f) over %f s\n", ...
%             x0, y0, th0, v0, x1, y1, th1, v1, T);

    d = [x0;
         y0;
         0];

    c = [f0 * cos(th0);
         f0 * sin(th0);
         v0 / f0];

    b = [(((3 * (x1 - x0)) / T) - 2 * f0 * cos(th0) - f1 * cos(th1)) / T;
         (((3 * (y1 - y0)) / T) - 2 * f0 * sin(th0) - f1 * sin(th1)) / T;
         (3 - ((2 * v0) / f0) - (v1 / f1)) / T];

    a = [(f0 * cos(th0) + f1 * cos(th1) - ((2 * (x1 - x0)) / T)) / (T^2);
         (f0 * sin(th0) + f1 * sin(th1) - ((2 * (y1 - y0)) / T)) / (T^2);
         ((v1 / f1) + (v0 / f0) - 2) / (T^2)];
end
