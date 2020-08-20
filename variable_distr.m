function [base_weights, pdf, cdf, variation] = variable_distr(diameter_sample, num_variations, plot_flg, three_means_only)


% Show probability plot only for 3 difference cases (left middle and right of the nerve)

% Variation of Guassian Wights
if three_means_only
    variation_range = [-.38 0.13 0.89];
else
    variation_range = linspace(-1, 1, num_variations);
end

[base_weights, pdf, cdf] = pdf_pan_paper(false);

leg = {};

for k = 1:length(variation_range)
    small_weight_var = base_weights.small * (1 + variation_range(k));
    med_wight_var = base_weights.medium * (1 - variation_range(k)*0.75);
    large_weight_var = base_weights.large * (1 - variation_range(k)*0.75);
    
    probability_sample = pdf(diameter_sample, small_weight_var, med_wight_var, large_weight_var);
    probability_sample = probability_sample./sum(probability_sample);
    
    variation.diameter_mean(k) = sum(diameter_sample.*probability_sample); % Mean diameter um
    variation.diameter_std(k) = std(diameter_sample, probability_sample); % st deviation
    variation.small_weight(k) = small_weight_var;
    variation.mediu_weight(k) = med_wight_var;
    variation.large_weight(k) = large_weight_var;
    
    if plot_flg && three_means_only
        plot(diameter_sample, probability_sample), hold on
        leg(end+1) = {['small_weight = ' num2str(small_weight_var) ', mean = ' num2str(variation.diameter_mean(k))]};
    end
end

if plot_flg && three_means_only
    legend(leg)
    title('Probability')
end

%% Mean and std over Variation of Guassian Wights

if plot_flg
    figure,
    subplot(211)
    plot(variation_range, variation.diameter_std), title('diameter std')
    subplot(212)
    plot(variation_range, variation.diameter_mean), title('diameter mean')
end