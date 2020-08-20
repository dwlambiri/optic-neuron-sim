function M = generate_model(varargin)
%% Parse inputs

f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
    function out = read_pair(name,def)
        ind = find(f(name));
        if ind, out = varargin{ind + 1};
        else out = def;
        end
    end

neuron_scale = 30;

nerve_r = read_pair('nerve_r', 1500);
bundle_r_range = read_pair('bundle_r_range', [250 300]);
min_bundles_dis = read_pair('min_bundles_dis', min(bundle_r_range)/10);
neuron_r_range = read_pair('neuron_r_range', [.5 3]*neuron_scale); % should be [.15 8] um
neuron_dens = read_pair('neuron_dens', 1);

file = read_pair('file', []);

%% Init

M = struct;

if ~isempty(file)
    file_full_addr = ['models\' file '.mat'];
    if exist(file_full_addr,'file') && ~any(f('rewrite'))
        fprintf('Reading existing model: %s\n', file);
        load(file_full_addr, 'M');
        if any(f('plot_model')), M.plot.model(); end
        return;
    end
end

%% Bundle GUI

% Create bundles
bund_g = [];
bundle_dens = 1;

function Update(~,~)
    t1 = str2num(h_nerv_r.String); %#ok<*ST2NM>
    t2 = str2num(h_min_dis.String);
    t3 = str2num(h_bundle_r_range.String);
    if isempty(t1) || ~all(size(t1) == [1 1])
        h_feedback.String = 'Bad Nerve Raduis!';
        return;
    end
    if isempty(t2) || ~all(size(t2) == [1 1])
        h_feedback.String = 'Bad Bundle Density!';
        return;
    end
    if isempty(t3) || ~all(size(t3) == [1 2])
        h_feedback.String = 'Bad Bundle Raduis Range!';
        return;
    end
    nerve_r = t1;
    min_bundles_dis = t2;
    bundle_r_range = t3;
    
    n_tries = 20;
    bund_g = [];
    while isempty(bund_g) && n_tries > 0
        bund_g = fill_circles(0, 'main_r', nerve_r, 'dens', bundle_dens, 'r_range', bundle_r_range, 'min_dis', min_bundles_dis);
        n_tries = n_tries - 1;
    end
    if ~n_tries
        h_feedback.String = 'Could not fit any bundles within the nerve, please adjust the parameters';
        return;
    end
    
    M.bund = bund_g;
    M.nerve_r = nerve_r;
    axis(ax_ui); cla;
    plot_model('no_neuron');
    h_feedback.String = sprintf('Model updated! (%d bundle(s))', size(bund_g, 2));
end

function h = pair_text(txt, pos, def)
    uicontrol('style','edit', 'Units','Normalized', 'position', [pos .08 .05], 'String', txt, 'Enable', 'off', 'BackgroundColor', [0.9400 0.9400 0.9400]);
    h = uicontrol('style','edit', 'Units','Normalized', 'position', [pos+[.08 0] .08 .05], 'String', num2str(def));
end
ok = 0;
    function OK(~,~)
        delete(fig_ui);
        ok = 1;
        drawnow;
    end

if any(f('GUI'))
    fig_ui = figure('units', 'normalized', 'toolbar','figure', 'Name', 'Bundle Settings', 'outerposition', [0.2 0.1 .6 .8]);
    ax_ui = axes('units', 'normalized', 'position', [0.1 0.1 .7 .8]);
    
    h_feedback = pair_text('Status: ', [.05 .92], []); h_feedback.Position = [.13 .92 .8 .05]; h_feedback.Enable = 'off';
    h_nerv_r = pair_text('Nerve Raduis:', [.82 .8], nerve_r);
    h_min_dis = pair_text('Min Clearance', [.82 .7], min_bundles_dis);
    h_bundle_r_range = pair_text('Bundle Raduis Range:', [.82 .6], bundle_r_range);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .5 .12 .05], 'String', 'Update', 'Callback', @Update);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .4 .12 .05], 'String', 'Done/Continue', 'Callback', @OK);
    Update(0,0);
    
    waitfor(fig_ui);
    if ~ok, return; end
else
    n_tries = 20;
    while isempty(bund_g) && n_tries > 0
        bund_g = fill_circles(0, 'main_r', nerve_r, 'dens', bundle_dens, 'r_range', bundle_r_range, 'min_dis', min_bundles_dis);
        n_tries = n_tries - 1;
    end
    if ~n_tries
        warning('Could not fit any bundles within the nerve, please adjust the parameters');
        return;
    end
end

%% Mesh Generation
fprintf('Generating model:\n');

% whole_geom = [[0 0 nerve_r]' bund_g];
n_bunds = size(bund_g, 2);
fprintf('\tNumber of bundles: %d\n', n_bunds);

neuron_g = cell(n_bunds,1);
expected_r_avg = zeros(n_bunds,1);

fprintf('\tNumber of neurons: ');
for k = 1:n_bunds
    expected_r_avg(k) = raduis_avg(bund_g(1:2,k)./nerve_r);
    neuron_g{k} = fill_circles(expected_r_avg(k), 'main_r', bund_g(3,k), 'dens', neuron_dens, 'r_range', neuron_r_range, 'min_dis', min(neuron_r_range)/2);
    fprintf('(%d) ', size(neuron_g{k}, 2));
    neuron_g{k}(1,:) = neuron_g{k}(1,:) + bund_g(1,k);
    neuron_g{k}(2,:) = neuron_g{k}(2,:) + bund_g(2,k);
end
fprintf('\n');

    function M = create_csg(M, varargin)
        f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
        if isfield(M, 'csg') && ~any(f('rewrite')), return, end
        fprintf('\tCreating CSG model... ');
        whole_geom = [[0 0 M.nerve_r]' M.bund M.neuron{:}];
        [M.csg.dl, M.csg.bt] = decsg([ones(1,size(whole_geom,2)); whole_geom]);
        fprintf('DONE\n');
    end

M.nerve_r = nerve_r;
M.neuron_r_range = neuron_r_range;
M.bund = bund_g;
M.neuron = neuron_g;
M.expected_r_avg = expected_r_avg;
M.create_csg = @create_csg;
M.plot.model = @plot_model;
M.plot.histogram = @plot_histogram;
M.plot.avg_r = @plot_avg_r;

M.file_full_addr = file_full_addr;
M.save = @(M) save(M.file_full_addr, 'M');

if any(f('plot_model')), plot_model(); end

    % finds the average neuron size based on statistics in: Pan et. al [2012]
    function r = raduis_avg(norm_cent)
        if ~isfield(M, 'r_dist')
            % Radius average in different regions
            v = [1.04 1.19	1.15    1.21    1.28    1.27
                .97  .95    1.01    1.16    1.24    1.26
                .93  .92    1.07    1.16    1.22    1.28
                .94  .93    1.08    1.19    1.18    1.17
                .92  1.02   1.02    1.14    1.16    1.24
                1.01 1.1    1.1     1.22    1.22    1.21  ] * neuron_scale;
            % imshow(-v, [], 'initialmagnification', 6000)
            stp = .01; % Interpolation step
            
            [x, y] = meshgrid(-1:(2/5):1);
            [xq, yq] = meshgrid(-1:stp:1); % Interpolate data
            r_avg = interp2(x,y,v,xq,yq);
            r_avg = imfilter(r_avg, fspecial('disk', .2/stp), 'symmetric');
            r_avg((xq.^2+yq.^2) > 1) = nan;
            M.r_dist.xq = xq;
            M.r_dist.yq = yq;
            M.r_dist.r_avg = r_avg;
        end
        
        dis2 = ((M.r_dist.xq-norm_cent(1)).^2 + (M.r_dist.yq-norm_cent(2)).^2);
        [~, ix] = min(dis2(:));
        
        r = M.r_dist.r_avg(ix);
    end

    function plot_avg_r
%         imagesc(M.r_dist.r_avg);
        surf(M.r_dist.xq, M.r_dist.yq, M.r_dist.r_avg, 'edgecolor', 'none'); view(0, 90);
        axis equal
    end


    function plot_model(varargin)
        ff = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
        no_neuron = ff('no_neuron'); varargin(no_neuron) = [];
        axis equal, hold on
        draw_circle([0 0], M.nerve_r, 200, varargin{:});
        draw_circle(M.bund(1:2,:)', M.bund(3,:)', 100, varargin{:});
        if ~any(no_neuron)
            for kk = 1:length(M.neuron)
                draw_circle(M.neuron{kk}(1:2,:)', M.neuron{kk}(3,:)', 20, varargin{:});
            end
        end
        drawnow
    end

    function plot_histogram(varargin)
        opts = {'Normalization', 'probability'};
        if isempty(varargin)
            fig = figure();
            title('Histogram: actual mean (expected mean)');
            plot_model('no_neuron');
            par_ax = gca;
            xl = get(gca, 'xlim');
            yl = get(gca, 'ylim');
            ax = zeros(length(M.neuron), 1);
            p = zeros(length(M.neuron), 4);
            for idx = 1:length(M.neuron)
                p(idx,:) = [M.bund(1:2,idx)'-M.bund(3,idx)*.7, M.bund(3,idx)*[1.4 1.4]];
                ax(idx) = axes();
                histogram(M.neuron{idx}(3,:), max(5, round(size(M.neuron{idx}, 2)/35)), opts{:});
                xlim(M.neuron_r_range);
                % set(ax(idx), 'XTickLabel', '');
                set(ax(idx), 'YTickLabel', '');
                text(.3, 1.1, sprintf('%.1f(%.1f)', mean(M.neuron{idx}(3,:)), M.expected_r_avg(idx)), 'Units', 'Normalized', 'Color', [.8 .2 .3]);
            end
            set(fig, 'SizeChangedFcn', @size_changed);
            size_changed(0, 0);
        else
            histogram(M.neuron{varargin{1}}(3,:), opts{:});
            xlim(M.neuron_r_range);
        end
        drawnow
        
        function size_changed(~,~)
            pos = utils.plotboxpos(par_ax);
            for id = 1:length(M.neuron)
                set(ax(id), 'position', [pos(1) + (p(id,1)-min(xl))/diff(xl)*pos(3),...
                    pos(2) + (p(id,2)-min(yl))/diff(yl)*pos(4),...
                    p(id,3)/diff(xl)*pos(3), p(id,4)/diff(yl)*pos(4)]);
            end
        end
    end

end

% Output: first row center_x
%    second row: center_y
%    third row: raduis
function geom = fill_circles(expected_r_avg, varargin)

Rand = @(siz, m) min(m) + rand(siz)*diff(m);

%% Input parsing
f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
    function out = read_pair(name,def)
        ind = find(f(name));
        if ind, out = varargin{ind + 1};
        else out = def;
        end
    end

main_r = read_pair('main_r', 20);
dens = read_pair('dens', 1);
r_range = read_pair('r_range', [1 5]);
min_dis = read_pair('min_dis', min(r_range));

n_max_circles = ceil((main_r / mean(r_range))^2); % normalize n_iters based on dens(ity)

%% Radius distribution
if expected_r_avg    
    n_max_circles = ceil((main_r / expected_r_avg)^2);
    mean_r_normalized = expected_r_avg / max(r_range);
    raduis_set = min(r_range) + skewed_distr(n_max_circles*2, mean_r_normalized) * diff(r_range);
    raduis_set = raduis_set(1:n_max_circles);
    raduis_set = sort(raduis_set, 'descend');
else
    raduis_set = Rand(n_max_circles, r_range);
end

%% Main loop
center = []; radius = [];
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
    r = raduis_set(k);
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

% creates skewed distribution of raduis of neurons
% siz: number of radii
% mean_r: mean radius
% x: array of size of siz with mean of mean_r and ~ lognormal distribution
function x = skewed_distr(siz, mean_r)

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

function draw_circle(center, raduis, res, varargin)

x = center(:,1);
y = center(:,2);

ang = (linspace(0,2*pi,res)); % min(round(.5*raduis), 15)
xp = raduis*cos(ang);
yp = raduis*sin(ang);
h = plot((x*ones(size(ang))+xp)',(y*ones(size(ang))+yp)', varargin{:});
for k = 1:length(h), h(k).ZData = ones(size(ang)); end

end