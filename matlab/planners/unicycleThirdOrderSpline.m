function [a, b, c, d] = unicycleThirdOrderSpline(x0, y0, th0, v0, x1, y1, th1, v1, T)
  % fprintf("Computing spline from (%f, %f, %f, %f) to (%f, %f, %f, %f)\n", ...
  %         x0, y0, th0, v0, x1, y1, th1, v1);

  d = [x0;
       y0];

  c = [v0 * cos(th0);
       v0 * sin(th0)];
  b = [(3 * (x1 - x0) - 3 * v0 * T * cos(th0) - T * (v1 * cos(th1) - v0 * cos(th0))) / T^2;
       (3 * (y1 - y0) - 3 * v0 * T * sin(th0) - T * (v1 * sin(th1) - v0 * sin(th0))) / T^2];
  a = [(v1 * cos(th1) - v0 * cos(th0) - 2 * b(1) * T) / (3 * T^2);
       (v1 * sin(th1) - v0 * sin(th0) - 2 * b(2) * T) / (3 * T^2)];
end
