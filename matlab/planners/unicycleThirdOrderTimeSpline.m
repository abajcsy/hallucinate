function [a, b, c, d] = unicycleThirdOrderTimeSpline(x0, y0, th0, v0, f0, x1, ...
                                                     y1, th1, v1, f1, T)
    % fprintf("Computing spline from (%f, %f, %f, %f) to (%f, %f, %f, %f)\n", ...
    %         x0, y0, th0, v0, x1, y1, th1, v1);

    % This spline has the following parametric form:
    %  x(p) = a(1) * p^3 + b(1) * p^2 + c(1) * p + d(1)
    %  y(p) = a(2) * p^3 + b(2) * p^2 + c(2) * p + d(2)
    %  p(t) = a(3) * t^3 + b(3) * t^2 + c(3) * t + d(3)

    d = [x0;
         y0;
         0];

    c = [v0 * cos(th0) / f0;
         v0 * sin(th0) / f0;
         f0];

    b = [(-(v1 * cos(th1)) / f1 + (3 * (x1 - x0)) / T - 2 * c(1)) / T;
         (-(v1 * sin(th1)) / f1 + (3 * (y1 - y0)) / T - 2 * c(2)) / T;
         (-f1 + 3 - 2 * c(3)) / T];

    a = [((v1 * cos(th1)) / f1 - c(1) - 2 * b(1) * T) / (3 * T^2);
         ((v1 * sin(th1)) / f1 - c(2) - 2 * b(2) * T) / (3 * T^2);
         (1 - c(3) - b(3) * T) / T^2];
end
