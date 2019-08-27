T = 1;
x0 = 0;
y0 = 0;
th0 = 0;
v0 = 0.1;
x1 = 1;
y1 = 1;
th1 = 0;
v1 = 0.1;

f0 = v0 + norm([x1 - x0; y1 - y0]);
f1 = v0 + norm([x1 - x0; y1 - y0]);
% f0 = 1;
% f1 = 1;
[a, b, c, d] = unicycleThirdOrderTimeSpline(x0, y0, th0, v0, f0, ...
                                            x1, y1, th1, v1, f1, ...
                                            T);

figure;
hold on;

pfunc = @(t) a(3) .* t.^3 + b(3) .* t.^2 + c(3) .* t + d(3);
xfunc = @(t) a(1) .* pfunc(t).^3 + b(1) .* pfunc(t).^2 + c(1) .* pfunc(t) + d(1);
yfunc = @(t) a(2) .* pfunc(t).^3 + b(2) .* pfunc(t).^2 + c(2) .* pfunc(t) + d(2);

fplot(xfunc, yfunc, [0 T]);

t = 0;
dt = 0.001;

x = getState(a, b, c, d, t);
xWithCtrl = getState(a, b, c, d, t)

xTri = drawTriangle([x(1); x(2)], x(3), 0.05);
xWithCtrlTri = drawTriangle([xWithCtrl(1); xWithCtrl(2)], xWithCtrl(3), 0.05);

while t <= T
    x = getState(a, b, c, d, t)
    ctrl = getControl(a, b, c, d, t);

    xWithCtrl = xWithCtrl + dt * [xWithCtrl(4) * cos(xWithCtrl(3)); ...
                        xWithCtrl(4) * sin(xWithCtrl(3)); ...
                        ctrl(1); ...
                        ctrl(2)]

    xTri = drawTriangle([x(1); x(2)], x(3), 0.05, xTri);
    xWithCtrlTri = drawTriangle([xWithCtrl(1); xWithCtrl(2)], xWithCtrl(3), ...
                                0.05, xWithCtrlTri);

    drawnow;
    pause(dt);
    t = t + dt;
end

hold off;

function state = getState(a, b, c, d, t)
    s = t;

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

function control = getControl(a, b, c, d, t)
    s = t;

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

function tri = drawTriangle(origin, rot, sideLength, lastTri)
    aLocal = [2*sideLength; 0; 1];
    bLocal = [0; -sideLength / sqrt(2); 1];
    cLocal = [0; sideLength / sqrt(2); 1];

    T = [cos(rot), -sin(rot), origin(1);
         sin(rot), cos(rot), origin(2);
         0, 0, 1];
    a = T * aLocal;
    b = T * bLocal;
    c = T * cLocal;

    if nargin > 3
        lastTri.set('XData', [a(1), b(1), c(1)]);
        lastTri.set('YData', [a(2), b(2), c(2)]);
        tri = lastTri;
    else
        tri = fill([a(1), b(1), c(1)], ...
                   [a(2), b(2), c(2)], ...
                   'b');
        set(tri, 'facealpha', 0.5);
    end
end
