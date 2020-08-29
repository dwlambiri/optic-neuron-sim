function x = skewed_distr(siz, median_r, mode_r, min_r, max_r)
% creates skewed distribution of radii of neurons
% siz: number of radii
% mean_r: mean radius
% x: array of size of siz with mean of mean_r and ~ lognormal distribution

% Initial setting of graph
%x_max = 100; % maximum x initial graph
mu = log(median_r-min_r);
%sig = sqrt(2*(log(mean_r-min_r)-mu));
sig = sqrt(mu - log(mode_r-min_r));


% Create initial distribution
x = lognrnd(mu,sig, [siz 1]);
%x(x > x_max) = [];
% histogram(x, 50, 'Normalization', 'probability');

% Normalize distribution range and adjust its mean
%x = (x-min(x))/(x_max  - min(x));
%x(round(siz*(mean_r-min_r) / mean(x)):end) = [];
%x = x * (mean_r-min_r) / mean(x);
% histogram(x, 50, 'Normalization', 'probability'); xlim([0 1]); Line(mean(x), 'V');
x = x+ min_r;
percentextra = 0.035;
p2p = ceil(percentextra*siz/100);
for i=40:60
    x = [ x ;i*ones(p2p,1)];
end
x(x > max_r) = [];
end
