function [base_weights, pdf, cdf] = pdf_pan_paper(print_flg)

%% Graph Data

epsil = 0.035;

diameter = [0;0.1;0.2;0.3;0.4;0.45;0.5;0.55;0.6;0.65;0.7;0.8;0.9;1;1.1;1.2;1.3;1.4;1.5;1.6;1.7;1.8;1.9;2;2.1;2.2;2.3;2.4;2.5;2.6;2.7;2.8;2.9;3;3.1;3.2;3.3;3.4;3.5;3.6;3.7;3.8;3.9;4;4.1;4.2;4.3;4.4;4.5;4.6;4.7;4.8;4.9;5;5.1;5.2;5.3;5.4;5.5;5.6;5.7;5.8;5.9];
probability = [0;0;0;0.7;1.25;2.95;4;5;6.5;8;9.4;9.6;8.7;7.7;6.7;5.7;4.85;4.1;3.65;3;2.65;2.25;2.1;1.85;1.4;1.3;1.2;1;0.9;0.8;0.65;0.55;0.55;0.45;0.4;0.35;0.3;0.25;0.25;0.2;0.1;0.1;0.1;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil];

probability = probability./sum(probability); % normalize probability

if print_flg
    fprintf('diameter mean = %f um\n', sum(diameter.*probability)); % 1.16 um
    fprintf('diameter std = %f um\n', std(diameter, probability)); % 0.70 um
end

%% Fitting Guassian 

% ============= Use "fit" for future contributions

diameter_sample = 0:0.001:6.1;

base_weights.small =        0.05674;%  (0.02847, 0.08502)
small_mean =                0.7301;%  (0.7167, 0.7435)
small_sig =                 0.2529;%  (0.212, 0.2938)
base_weights.medium =       0.03945;%  (0.01806, 0.06084)
med_mean =                  1.076;%  (0.888, 1.265)
med_sig =                   0.4416;%  (0.2527, 0.6304)
base_weights.large =        0.01741;%  (0.01206, 0.02277)
large_mean =                1.73;%  (1.367, 2.094)
large_sig =                 0.9244;%  (0.7075, 1.141)

pdf = @(x, w1, w2, w3)    w1*exp(-((x-small_mean)/small_sig).^2) + ...
                          w2*exp(-((x-med_mean)/med_sig).^2) + ...
                          w3*exp(-((x-large_mean)/large_sig).^2);

probability_sample = pdf(diameter_sample, base_weights.small, base_weights.medium, base_weights.large);
probability_sample = probability_sample./sum(probability_sample); % normalize probability

if print_flg
    fprintf('fit diameter mean = %f um\n', sum(diameter_sample.*probability_sample)); % 1.19 um
    fprintf('fit diameter std = %f um\n', std(diameter_sample, probability_sample)); % 0.59 um
end

cdf = @(x, w1, w2, w3) (w1 * normcdf(x, small_mean, small_sig) +...
                        w2 * normcdf(x, med_mean, med_sig) + ...
                        w3 * normcdf(x, large_mean, large_sig))...
                        /(w1 + w2 + w3);

