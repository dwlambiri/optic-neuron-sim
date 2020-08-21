function x = skewed_distr(siz, mean_r)
% creates skewed distribution of radii of neurons
% siz: number of radii
% mean_r: mean radius
% x: array of size of siz with mean of mean_r and ~ lognormal distribution

% Initial setting of garph
E = 5; % Expected value
x_max = 30; % maximum x initial graph
mu = 1.4;

% Create initial distribution
sig = sqrt(2*(log(E) - mu)); % E = exp(mu + sig^2/2);
x = lognrnd(mu,sig, [siz 1]);
x(x > x_max) = [];
% histogram(x, 50, 'Normalization', 'probability');

% Normalize distribution range and adjust its mean
x = (x-min(x))/(x_max - min(x));
x(round(siz*mean_r / mean(x)):end) = [];
x = x * mean_r / mean(x);
% histogram(x, 50, 'Normalization', 'probability'); xlim([0 1]); Line(mean(x), 'V');

end
