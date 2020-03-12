clc
clear all

%% Priors.
prior_g1 = 0.5;
prior_g2 = 0.5;

%% Means (and variance) of gaussians.
mu_g1 = -pi/4;
mu_g2 = pi/4;
var = pi/8;

%% Likelihood of different actions under gaussian.
down = -pi/2;
pdown_g1 = gaussian_obs(down, mu_g1, var);
pdown_g2 = gaussian_obs(down, mu_g2, var);

%pdown_g1 = normpdf(down, mu_g1, var)
%pdown_g2 = normpdf(down, mu_g2, var)

u_g1_opt = -pi/4;
pu_g1_opt_g1 = gaussian_obs(u_g1_opt, mu_g1, var);
pu_g1_opt_g2 = gaussian_obs(u_g1_opt, mu_g2, var);

%% Calculate posterior for u=DOWN.
pg1_down = (pdown_g1 * prior_g1) / ...
            (pdown_g1 * prior_g1 + pdown_g2 * prior_g2)
        
pg2_down = (pdown_g2 * prior_g2) / ...
            (pdown_g1 * prior_g1 + pdown_g2 * prior_g2)

%% Calculate posterior for u=OPT_G1_ACTION.
pg1_opt = (pu_g1_opt_g1 * prior_g1) / ...
            (pu_g1_opt_g1 * prior_g1 + pu_g1_opt_g2 * prior_g2)
        
pg2_opt = (pu_g1_opt_g2 * prior_g2) / ...
            (pu_g1_opt_g1 * prior_g1 + pu_g1_opt_g2 * prior_g2)

function pu = gaussian_obs(u, mu, var)
    % Normalizer.
    z = 1/(sqrt(2 * pi * var)); 
    pu = z * exp(-(0.5/var) * (u - mu)^2);
end