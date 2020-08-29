
close all

study_three_means_only = false;

num_variations = 1000;
num_diameter_samples = 1000;
num_cdf_samples = 70;

min_abs_dia = 0.19;
max_abs_dia = 6.87;

min_mean_dia = 0.92;
max_mean_dia = 1.28;

diameter_sample = linspace(min_abs_dia, max_abs_dia, num_diameter_samples);

[base_weights, pdf, cdf, variation] = variable_distr(diameter_sample, num_variations, false, study_three_means_only);

out_of_bound = variation.diameter_mean < min_mean_dia | variation.diameter_mean > max_mean_dia;

variation.small_weight(out_of_bound) = [];
variation.mediu_weight(out_of_bound) = [];
variation.large_weight(out_of_bound) = [];
variation.diameter_mean(out_of_bound) = [];

cdf_interpolated = linspace(0,1,num_cdf_samples);

n_dias = length(variation.diameter_mean);

cdf_tbl = single(zeros(n_dias, num_cdf_samples+1));

cdf_tbl(:,1) = variation.diameter_mean; % First column: average diameter

possible_underflow_locaions = 1:num_cdf_samples < num_cdf_samples * .1;
possible_overrflow_locaions = 1:num_cdf_samples > num_cdf_samples * .9;

for k = 1:n_dias
    
    cdf_k = cdf(diameter_sample,...
            variation.small_weight(k),...
            variation.mediu_weight(k),...
            variation.large_weight(k));
            
    diameter_interpolated = interp1(cdf_k, diameter_sample, cdf_interpolated);
    
    out_of_bounds = isnan(diameter_interpolated) > 0;
    diameter_interpolated(out_of_bounds & possible_underflow_locaions) = min_abs_dia;
    diameter_interpolated(out_of_bounds & possible_overrflow_locaions) = max_abs_dia;
    
    cdf_tbl(k,2:end) = diameter_interpolated./2; % >>>>> RADIUS
    
    if study_three_means_only || k == 1 || k == n_dias/2 || k == n_dias
        figure
        %diameter_interpolated
        %cdf_interpolated
        plot(diameter_interpolated, cdf_interpolated)
        title(sprintf('CDF for mean = %.2f um', variation.diameter_mean(k)))
        xlabel('Axon diameter(um)');
        ylabel('Cumulative Probability');       
    
        figure
        a = diff(cdf_interpolated)./diff(diameter_interpolated);
        diameter_interpolated(length(diameter_interpolated)) = [];
        b = diameter_interpolated;
        plot(b, a)
        title(sprintf('PDF for mean = %.2f um', variation.diameter_mean(k)))
        xlabel('Axon diameter(um)');
        ylabel('Probability');    
    end
    
end

%% Write the table to file
if ~study_three_means_only
%     whos cdf_tbl % sse size of data in bytes
%     fid = fopen('axon_radius_cdf.tbl', 'w');
%     fwrite(fid, size(cdf_tbl), 'uint32');
%     fwrite(fid, cdf_tbl, 'single');
%     fclose(fid);
end
