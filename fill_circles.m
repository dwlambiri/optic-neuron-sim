% Output: first row center_x
%    second row: center_y
%    third row: raduis
function geom = fill_circles(expected_r_avg, main_r, dens, r_range, min_dis)

Rand = @(siz, m) min(m) + rand(siz)*diff(m);
del = @(x) fprintf(repmat('\b',1,x));


%% Radius distribution
if expected_r_avg    
    n_max_circles = ceil((main_r / expected_r_avg)^2);
    %mean_r_normalized = expected_r_avg / max(r_range);
    radius_set = skewed_distr(n_max_circles, expected_r_avg, min(r_range), max(r_range));
    n_max_circles = size(radius_set,1);
    radius_set = radius_set(1:n_max_circles);
    radius_set = sort(radius_set, 'descend');
else
    n_max_circles = ceil((main_r / mean(r_range))^2); % normalize n_iters based on dens(ity)
    radius_set = Rand(n_max_circles, r_range);
end

%% Main loop
if expected_r_avg  
    center = [0 0]; 
    radius(1) = main_r/10;
    sumrad2 = (main_r/20)^2;
else
    center = [];
    radius = [];
    sumrad2= 0;
end
    
nDispChar = 0; last_disp_num = 0;
tic
for k = 1:n_max_circles
    % display
    new_disp_num = round(100*k/n_max_circles);
    if new_disp_num ~= last_disp_num
        del(nDispChar); nDispChar = fprintf('%d of %d\t', k, n_max_circles);
        last_disp_num = new_disp_num;
    end

    %% Create a random circle
    r = radius_set(k);

    if ~isempty(radius)
        if (main_r^2 - sumrad2 - r^2)<=0 
            continue;
        end
    end
    emptySpace = main_r^2 - sumrad2;
    if emptySpace == 0
        break;
    end
    n_tries = ceil((main_r^2/emptySpace)^3);
    %n_tries = (dens + k);
    center_set = Rand([n_tries 2], (main_r-r-min_dis)*[-1 1]);
    for kk = 1:n_tries % tring to fit the circle with raduis r
        c = center_set(kk,:);
        % Check if non overlapping
        if (norm(c) + r + min_dis < main_r)
            if isempty(radius)
                %tic
                center = [center; c];
                radius = [radius; r];
                sumrad2 = sumrad2 + r^2;
                %toc
                break;
            else
                C2Cdis = sqrt((center(:,1) - c(1)).^2 + (center(:,2) - c(2)).^2);
                if all(C2Cdis > radius + r + min_dis)
                    %tic
                    center = [center; c];
                    radius = [radius; r];
                    sumrad2 = sumrad2 + r^2;
                    %toc
                    break;
                end
            end
        end
    end
end
toc
del(nDispChar);
if expected_r_avg  
    if length(center(:,1)) >= 1
        %sumrad2 = sumrad2 - radius(1)^2;
        center(1,:) = [];
        radius(1) = [];
    end
end
geom = [center radius]';
fprintf("Ratio of filled space %0.2f\n",sumrad2/main_r^2);
fprintf('Expected radius %f min rad %f max rad %f\n', expected_r_avg, min(r_range), max(r_range));
fprintf("Average radius %f, min %f, max %f std %f\n", mean(radius), min(radius), max(radius), std(radius));
end

