% Output: first row center_x
%    second row: center_y
%    third row: raduis
function geom = fill_circles(expected_r_avg, main_r, dens, r_range, min_dis, hard, obstacles, circleCenters)

Rand = @(siz, m) min(m) + rand(siz)*diff(m);
del = @(x) fprintf(repmat('\b',1,x));

if ~isempty(circleCenters)
    size(circleCenters) 
end
%% Radius distribution
if expected_r_avg    
    n_max_circles = ceil((main_r / expected_r_avg)^2);
    %mean_r_normalized = expected_r_avg / max(r_range);
    radius_set = skewed_distr(n_max_circles, expected_r_avg, min(r_range), max(r_range));
    n_max_circles = size(radius_set,1);
    radius_set = radius_set(1:n_max_circles);
    radius_set = sort(radius_set, 'descend');
else
    %n_max_circles = ceil((main_r / mean(r_range))^2); % normalize n_iters based on dens(ity)
    %radius_set = Rand(n_max_circles, r_range);
    v1 = mean(r_range);
    n_max_circles = ceil((main_r / v1)^2);
    %mean_r_normalized = expected_r_avg / max(r_range);
    radius_set = skewed_distr(n_max_circles, v1, min(r_range), max(r_range));
    n_max_circles = size(radius_set,1);
    radius_set = radius_set(1:n_max_circles);
    radius_set = sort(radius_set, 'descend');
end

    function [n, centers, radii] = make_capillary(line, r_small)
        x = line(1);
        y = line(2);
        x1 = line(3);
        y1 = line(4);
        n = floor(sqrt((x-x1)^2+(y-y1)^2)/(2*r_small));
        if n < 1
            centers = [];
            radii = [];
            return;
        end
        fprintf("r+small = %f N= %d\n", r_small, n);
        centers = [linspace(x, x1, n)' linspace(y, y1, n)']; 
        radii = linspace(r_small, r_small, n)';
        
    end
rem = 0;
%% Main loop
sumrad2 = 0;
if obstacles  
    ratio_1 = 6.2143;
    ratio_2 = 7.25;
    

    center = [(0+main_r/ratio_1) 0]; % create artery
    radius = main_r/ratio_1;
    %sumrad2 = (main_r/ratio_1)^2;
    
    center = [center; (0-main_r/ratio_2) 0]; % create vein
    radius = [radius; main_r/ratio_2];
    %sumrad2 = sumrad2 + (main_r/ratio_2)^2;
    
    m=400; 
    [vv, cc] = voronoin([2*main_r*(rand(1,m)-0.5)' 2*main_r*(rand(1,m)-0.5)']);
    
    total = 0;
    for q=1:length(cc)
        total = total + length(cc{q})-1;
    end
    
    lines = zeros(4,total);
    good = vv < main_r & vv > -main_r;
    good = good(:,1) & good(:,2);
    iter = 1;
    for q=1:length(cc)
        face = cc{q};
        for w=1:length(face)-1
           if ~good(face(w)) || ~good(face(w+1))
               continue;
           end 
           lines(1,iter) = vv(face(w),1);
           lines(2,iter) = vv(face(w),2);
           lines(3,iter) = vv(face(w+1),1);
           lines(4,iter) = vv(face(w+1),2);
           iter = iter+1;
        end
        if good(face(1)) && good(face(length(face)))
           lines(1,iter) = vv(face(length(face)),1);
           lines(2,iter) = vv(face(length(face)),2);
           lines(3,iter) = vv(face(1),1);
           lines(4,iter) = vv(face(1),2);
           iter = iter+1;
        end
    end

    rem = 2;
    
    num = iter;
    obstacle_r = (50-m/20)*rand(1,num);
    
    fprintf("Number of obstacles = %d\n", num);
    
    for i =1:num-1
        
        if obstacle_r(i) < 5
            continue;
        end
        [n, b, g] = make_capillary(lines(:,i), obstacle_r(i));
        if isempty(b) || isempty(g)
            continue;
        end
        for k =1:n-1
                if isempty(center)
                    center = b(1,:);
                    radius = g(1);
                else
                    %fprintf("i = %d, k=%d\n",i, k);
                    C2Cdis = sqrt((center(:,1) - b(k,1)).^2 + (center(:,2) - b(k,2)).^2);
                    if all(C2Cdis > radius + g(k))
                        %tic
                        center = [center; b(k,:)];
                        radius = [radius; g(k)];
                        rem = rem+1;
                        %toc
                    end
              end
        end
    end
    
else
    center = [];
    radius = [];
    sumrad2= 0;
end

fprintf(" REM = %d\n", rem);
    
nDispChar = 0; last_disp_num = 0;
if hard > 4 
   hard = 4;
end
if hard < 0
   hard = 0;
end
hard = mod(ceil(hard),4)+1;

zones = [-1 -2/3 -1/3 0 1/3 2/3];

tic
countNeurons = 0;
%maxNeurons = ceil(160000* (main_r/10000)^2);
maxNeurons = 10^6;
stopPlacing = false;
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

    n_tries = ceil((main_r^2/emptySpace)^hard);
    %fprintf(" Hard = %d Tries = %d\n", hard, n_tries);
    %n_tries = (dens + k);
    center_set = Rand([n_tries 2], (main_r-r)*[-1 1]);
    for kk = 1:n_tries % tring to fit the circle with raduis r
        c = center_set(kk,:);
        % Check if non overlapping
        if (norm(c) + r  < main_r)
            if isempty(radius)
                if ~isempty(circleCenters)
                    xt = floor(main_r+c(1));
                    yt = floor(main_r-c(2));
                    if xt < 1
                        xt = 1;
                    end
                    if yt < 1
                        yt = 1;
                    end
                    if circleCenters(yt,xt) == 1
                        %fprintf('hit empty\n');
                        continue;
                    end
                end
                %tic
                center = [center; c];
                radius = [radius; r];
                sumrad2 = sumrad2 + r^2;
                countNeurons = countNeurons + 1;
                %toc
                break;
            else
                C2Cdis = sqrt((center(:,1) - c(1)).^2 + (center(:,2) - c(2)).^2);
                if ~isempty(circleCenters)
                    xt = floor(main_r+c(1));
                    yt = floor(main_r-c(2));
                    if xt < 1
                        xt = 1;
                    end
                    if yt < 1
                        yt = 1;
                    end
                    if circleCenters(yt,xt) == 1
                        %fprintf('hit empty\n');
                        continue;
                    end
                end
                if all(C2Cdis > radius + r + min_dis*rand())
                    %tic
                    center = [center; c];
                    radius = [radius; r];
                    sumrad2 = sumrad2 + r^2;
                    countNeurons = countNeurons + 1;
                    if countNeurons > maxNeurons
                        stopPlacing = true;
                    end
                    %toc
                    break;
                end 
            end
        end
    end
    if stopPlacing == true
         break;
    end
end
toc
del(nDispChar);
ra = 0;
if obstacles 
    if length(center(:,1)) >= 1
        ra = radius(2)^2 + radius(1)^2;
        %sumrad2 = sumrad2 - ra;
        
        for i = 1:rem
           center(1,:) = [];
           radius(1) = []; 
        end


    end
end
if isempty(center) 
    c = [0 0];
    r = min(r_range);
    center = [center; c];
    radius = [radius; r];
    sumrad2 = sumrad2 + r^2;
end
geom = [center radius]';
fprintf("Ratio of filled space %0.2f\n",sumrad2/(main_r^2-ra));
fprintf('Expected radius %f min rad %f max rad %f\n', expected_r_avg, min(r_range), max(r_range));
fprintf("Average radius %f, min %f, max %f std %f\n", mean(radius), min(radius), max(radius), std(radius));
end

