
goals = {[2,2], [2,-2]};
us = [pi/3, 0, -pi/3, -(2*pi)/3, -pi, (2*pi)/3];
x0 = [0, 0];
sigma = pi/4;
goalIdx = 1;

sz = 1;
total_p = zeros(sz,sz);
for u=us
    uarr = u * ones(sz,sz);
    xs = zeros(sz,sz);
    ys = zeros(sz,sz);
    u
    prob = Pu_given_x_g_normalized(uarr, xs, ys, goals, goalIdx, sigma, us)
    total_p = total_p + prob;
end
total_p

function prob = Pu_given_x_g_normalized(u, x, y, goals, goalIdx, sigma, us)
    % Get the lower and upper bounds to integrate the Gaussian 
    % PDF over.
    uopt = atan2(goals{goalIdx}(2)- y, goals{goalIdx}(1) - x);

    % minimum angular distance between current control (u) and uopt
    diff = abs(angdiff(u, uopt));

    pd1 = makedist('Normal','mu',0,'sigma', sigma);
    truncpd = truncate(pd1, -pi, pi);

    % Get all the probabilities
    unnormalized_probs = pdf(truncpd, diff);

    % normalize.
    norm = zeros(size(u));
    for otheru=us
        otheru_arr = otheru * ones(size(u));
        otherdiff = abs(angdiff(otheru_arr, uopt));
        new_prob = pdf(truncpd, otherdiff);
        norm = norm + new_prob;
    end

    prob = unnormalized_probs ./ norm;
end