% Initial state.
x0 = 0.;
y0 = 0.;
th0 = 0.;
v0 = 0.1;

% Final state.
% x1 = 1.;
x1 = 0.;
y1 = 1.;
th1 = 0.;
% th1 = pi/2.;
% th1 = pi/1.5;
% th1 = 3*pi/2.;
% th1 = pi;
v1 = 0.1;
% v1 = 0.5;
% v1 = 1.5;

T = 1.;
% T = 5.;

% Compute the spline.
[a, b, c, d] = unicycleThirdOrderSpline(x0, y0, th0, v0, x1, y1, th1, v1, T);

xfunc = @(t) a(1) .* t.^3 + b(1) .* t.^2 + c(1) .* t + d(1);
yfunc = @(t) a(2) .* t.^3 + b(2) .* t.^2 + c(2) .* t + d(2);

vxfunc = @(t) 3 .* a(1) .* t.^2 + 2 .* b(1) .* t + c(1);
vyfunc = @(t) 3 .* a(2) .* t.^2 + 2 .* b(2) .* t + c(2);

fprintf("x(T) = %f\n", xfunc(T));
fprintf("y(T) = %f\n", yfunc(T));
fprintf("vx(T) = %f\n", vxfunc(T));
fprintf("vy(T) = %f\n", vyfunc(T));

figure;
hold on;
fplot(xfunc, yfunc, [0, T]);
hold off;

