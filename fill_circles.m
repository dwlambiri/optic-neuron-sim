% Output: first row center_x
%    second row: center_y
%    third row: raduis
function geom = fill_circles(median_r_value, main_r, dens, r_range, min_dis, hard, obstacles, circleCenters, zoned)

Rand = @(siz, m) min(m) + rand(siz)*diff(m);
RandSq = @(siz, m) [ min(m(1,:))+rand(siz(1),1)*diff(m(1,:)) min(m(2,:))+rand(siz(1),1)*diff(m(2,:))];

max_n_tries = 50;

del = @(x) fprintf(repmat('\b',1,x));

if ~isempty(circleCenters)
    %size(circleCenters) 
    obstacles = false;
end

if zoned
    zoneSplit = [2/3 1; 1/3 2/3; 0 1/3; -1/3 0; -2/3 -1/3; -1 -2/3];
    yvalue = [sqrt(5/9); sqrt(8/9); 1 ; 1; sqrt(8/9); sqrt(5/9)];
    zspace = [1 0.8 0.7 0.7 1 1];
    nz = size(zoneSplit,1);
else
    zoneSplit = [-1 1];
    nz = 1;
    zspace = 1;
end

n_max_circles = zeros(1,nz);
%% Radius distribution
if median_r_value    
    %n_max_circles = ceil((main_r / mode_r_value)^2/6);
    % max circles = (114*1000*Area + 218000)*reduction
    %variability = 1 + (mod((rand()*100),20))/100;
    variability = 1.5;
    n_max = ceil(variability*(114000*3.1415*(main_r/1000)^2 + 218000)*(1/10)^2);
    if zoned
        per_zone_neurons = [.11 .18 .21 .21 .18 .11];
        divider = 5;
    else
        per_zone_neurons = 1;
        divider = 1;
    end
    gen_n = ceil(n_max/divider);
    radius_set = zeros(gen_n,nz);
    for i=1:nz
        %skewed_distr(siz, median_r, mode_r, min_r, max_r)
        if nz == 1
            radius_set_tmp = skewed_distr(gen_n, 2*median_r_value, 2*(median_r_value-1), 2*min(r_range), 2*(max(r_range)-3*(i-1)));
        else
            % use the formula from Pan figure 5 but for median not average!
            % median *diameter* - 0.0639*i + 0.8839 (um) -> for pixels multiply by 10
            % mode is minus 0.2 um is minux 2 pixels
            radius_set_tmp = skewed_distr(gen_n, 0.639*i +8.838, 0.639*i +6.838 , 2*min(r_range), 2*(max(r_range)-3*(i-1)));
        end
        radius_set_tmp = radius_set_tmp/2;
        radius_set_tmp = radius_set_tmp(1:length(radius_set_tmp));
        radius_set_tmp = sort(radius_set_tmp, 'descend');
        radius_set(1:length(radius_set_tmp),i) = radius_set_tmp;
        n_max_circles(i) = min(floor(n_max*per_zone_neurons(i)), length(radius_set_tmp));
    end
    %size( radius_set)
    %size( n_max_circles)
    %mean_r_normalized = expected_r_avg / max(r_range);
    %calculate distribution on diameters then divide by 2 for radius
else
    %n_max_circles = ceil((main_r / mean(r_range))^2); % normalize n_iters based on dens(ity)
    %radius_set = Rand(n_max_circles, r_range);
    v1 = mean(r_range);
    n_max_circles = ceil((main_r / v1)^2);
    %mean_r_normalized = expected_r_avg / max(r_range);
    radius_set = skewed_distr(n_max_circles, 2*v1,2*v1, 2*min(r_range), 2*max(r_range));
    radius_set = radius_set/2;
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
    for zoneIter=1:length(cc)
        total = total + length(cc{zoneIter})-1;
    end
    
    lines = zeros(4,total);
    good = vv < main_r & vv > -main_r;
    good = good(:,1) & good(:,2);
    iter = 1;
    for zoneIter=1:length(cc)
        face = cc{zoneIter};
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
        for circleIter =1:n-1
                if isempty(center)
                    center = b(1,:);
                    radius = g(1);
                else
                    %fprintf("i = %d, k=%d\n",i, k);
                    C2Cdis = sqrt((center(:,1) - b(circleIter,1)).^2 + (center(:,2) - b(circleIter,2)).^2);
                    if all(C2Cdis > radius + g(circleIter))
                        %tic
                        center = [center; b(circleIter,:)];
                        radius = [radius; g(circleIter)];
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

tic
countNeurons = 0;
%maxNeurons = ceil(160000* (main_r/10000)^2);
maxNeurons = 10^6;
stopPlacing = false;

sumradTemp = sumrad2;

for zoneIter = 1:nz
  zz = zoneSplit(zoneIter,:);
  sumrad2 = 0;
  
  for circleIter = 1:n_max_circles(zoneIter)
    % display
    new_disp_num = round(100*circleIter/n_max_circles(zoneIter));
    if new_disp_num ~= last_disp_num
        del(nDispChar); nDispChar = fprintf('%d of %d\t', circleIter, n_max_circles(zoneIter));
        last_disp_num = new_disp_num;
    end

    %% Create a random circle
    if nz == 1
        r = radius_set(circleIter);
    else
       r = radius_set(circleIter, zoneIter); 
    end

    if ~isempty(radius)
       if (main_r^2*zspace(zoneIter)/nz - sumrad2 - r^2)<=0 
           continue;
       end
    end
    
    emptySpace = main_r^2*zspace(zoneIter)/nz - sumrad2;
    if emptySpace == 0
        break;
    end

    if zoned == false
        n_tries = ceil((main_r^2*zspace(zoneIter)/nz/emptySpace)^(max(0,hard)));
    else
        n_tries = max_n_tries;
    end
    if n_tries > max_n_tries
        n_tries = max_n_tries;
    end
    
    %fprintf(" Hard = %d Tries = %d\n", hard, n_tries);
    %n_tries = (dens + k);
    if nz == 1
        center_set = Rand([n_tries 2], (main_r-r)*[-1 1]);
    else
        center_set = RandSq([n_tries 2], (main_r-r)*[zz; -1*yvalue(zoneIter) yvalue(zoneIter)]);
    end
    for tryIter = 1:n_tries % tring to fit the circle with radius r
        c = center_set(tryIter,:);
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
geom = [center radius zeros(length(radius),1)]';
fprintf("Ratio of filled space %0.2f\n",sumrad2/(main_r^2-ra));
fprintf('Expected Radius Median %f, Mode %f, Min %f, Max %f\n', median_r_value, median_r_value-1, min(r_range), max(r_range));
[N,P]= hist(radius, 60);
idx = find(N==max(N));
fprintf("Actual Radius   Median  %f, Mode %f, Min %f, Max %f std %f\n", median(radius), P(idx), min(radius), max(radius), std(radius));
fprintf('Total Placed Neurons %d\n', length(radius));
end

