function x = skewed_distr(siz, median_r, mode_r, min_r, max_r)
% creates skewed distribution of radii of neurons
% siz: number of radii
% median_r: median of distribution
% mode_r: mode of distribution
% min_r: min of distribution
% max_r: max of distribution


mu = log(median_r);
%sig = sqrt(2*(log(mean_r-min_r)-mu));
sig = sqrt(mu - log(mode_r));
% Create initial distribution
x = lognrnd(mu,sig, [siz 1]);

x = x+ (min_r - min(x));
percentextra = 0.035;

% add to the tail do the distribution
% as per Pan study the distribution has a long flat tail
% thus it is not a strict lognormal distribution
p2p = ceil(percentextra*siz/100);
for i=40:60
    x = [ x ;i*ones(p2p,1)];
end
x(x > max_r) = [];
end
