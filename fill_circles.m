% Output: first row center_x
%    second row: center_y
%    third row: raduis
function geom = fill_circles(expected_r_avg, main_r, dens, r_range, min_dis)

Rand = @(siz, m) min(m) + rand(siz)*diff(m);


n_max_circles = ceil((main_r / mean(r_range))^2); % normalize n_iters based on dens(ity)

%% Radius distribution
if expected_r_avg    
    n_max_circles = ceil((main_r / expected_r_avg)^2);
    mean_r_normalized = expected_r_avg / max(r_range);
    radius_set = min(r_range) + skewed_distr(n_max_circles*2, mean_r_normalized) * diff(r_range);
    radius_set = radius_set(1:n_max_circles);
    radius_set = sort(radius_set, 'descend');
else
    radius_set = Rand(n_max_circles, r_range);
end

%% Main loop
center = []; 
radius = [];

del = @(x) fprintf(repmat('\b',1,x));
nDispChar = 0; last_disp_num = 0;
for k = 1:n_max_circles
    % display
    new_disp_num = round(100*k/n_max_circles);
    if new_disp_num ~= last_disp_num
        del(nDispChar); nDispChar = fprintf('%.0f%%\t', new_disp_num);
        last_disp_num = new_disp_num;
    end

    %% Create a random circle
    r = radius_set(k);
    n_tries = (dens + k);
    center_set = Rand([n_tries 2], (main_r-r-min_dis)*[-1 1]);
    for kk = 1:n_tries % tring to fit the circle with raduis r
        c = center_set(kk,:);
        % Check if non overlapping
        if (norm(c) + r + min_dis < main_r)
            if isempty(radius)
                center = [center; c];
                radius = [radius; r];
                break;
            else
                C2Cdis = sqrt((center(:,1) - c(1)).^2 + (center(:,2) - c(2)).^2);
                if all(C2Cdis > radius + r + min_dis)
                    center = [center; c];
                    radius = [radius; r];
                    break;
                end
            end
        end
    end
end
del(nDispChar);

geom = [center radius]';

end

