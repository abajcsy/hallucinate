% Initial state.
x0 = 0.;
y0 = 0.;
th0 = 0.;
v0 = 0.1;

% Final state.
% x1 = 1.;
x1 = 0.;
y1 = 1.5;
% th1 = 0.;
th1 = pi/2.;
% th1 = pi/1.5;
% th1 = 3*pi/2.;
% th1 = pi;
v1 = 0.1;
% v1 = 0.5;
% v1 = 1.5;

T = 1.;
% T = 5.;

f0 = 1;
f1 = 1;

% Compute the spline.
[a, b, c, d] = unicycleThirdOrderTimeSpline(x0, y0, th0, v0, f0, x1, y1, th1, ...
                                            v1, f1, T);

pfunc = @(t) a(3) .* t.^3 + b(3) .* t.^2 + c(3) .* t + d(3);
vpfunc = @(t) 3 .* a(3) .* t.^2 + 2 .* b(3) .* t + c(3);

xfunc = @(t) a(1) .* pfunc(t).^3 + b(1) .* pfunc(t).^2 + c(1) .* pfunc(t) + d(1);
yfunc = @(t) a(2) .* pfunc(t).^3 + b(2) .* pfunc(t).^2 + c(2) .* pfunc(t) + d(2);

vxfunc = @(t) (3 .* a(1) .* pfunc(t).^2 + 2 .* b(1) .* pfunc(t) + c(1)) * vpfunc(t);
vyfunc = @(t) (3 .* a(2) .* pfunc(t).^2 + 2 .* b(2) .* pfunc(t) + c(2)) * vpfunc(t);

fprintf("x(T) = %f\n", xfunc(T));
fprintf("y(T) = %f\n", yfunc(T));
fprintf("vx(T) = %f\n", vxfunc(T));
fprintf("vy(T) = %f\n", vyfunc(T));

figure;
hold on;
fplot(xfunc, yfunc, [0, T]);

xsamples = [];
ysamples = [];

numSamples = 25;
Tsample = T / numSamples;
for t = 1:numSamples
    xsamp = xfunc(t * Tsample);
    ysamp = yfunc(t * Tsample);

    xsamples = [xsamples, xsamp];
    ysamples = [ysamples, ysamp];

    vxsamp = vxfunc(t * Tsample);
    vysamp = vyfunc(t * Tsample);
    vsamp = sqrt(vxsamp^2 + vysamp^2);
    fprintf('v(%f) = %f\n', t * Tsample, vsamp);
end

scatter(xsamples, ysamples);

hold off;
