function M = generate_model(varargin)
% GENERATE_MODEL
% parameters:
%   Type, Value
%       neuron_scale
%       nerve_r
%       bundle_r_range
%       min_bundle_dis
%       neuron_r_range
%       neuron_dens
%       file
%   Binary
%       rewrite
%       GUI
%       plot_model

M = struct;

%% Parse inputs

f = @(g) (cellfun(@(x) ischar(x) && strcmp(x,g), varargin));

    function out = read_pair(name,def)
        ind = find(f(name));
        if ind
            out = varargin{ind + 1};
        else
            out = def;
        end
    end

neuron_scale = read_pair('neuron_scale', 30);
nerve_r = read_pair('nerve_r', 1500);
bundle_r_range = read_pair('bundle_r_range', [1499 1500]);
min_bundles_dis = read_pair('min_bundles_dis', 2);
%neuron_r_range = read_pair('neuron_r_range', [.75 12]*neuron_scale); % should be [.15 7] um
neuron_r_range = [0.9 39]; % because of 10 pixels/um and radius
neuron_dens = read_pair('neuron_dens', 1);
refine = read_pair('refine', 0);

% PROPAGATION ALGO PARAMETERS
init_insult = read_pair('init_insult', [0 0]);
medium_speed = read_pair('medium_speed', .02);
bundle_speed = read_pair('bundle_speed', .02);
neuron_speed = read_pair('neuron_speed', 1);
all_edge_delay = read_pair('edge_delay', 2);

obstaclesOnBundles = true;


file = read_pair('file', []);

imageFile = read_pair('image', []);
oni = [];


if isempty(imageFile) 
    imageFile = 'optic-nerve.png';
end

if ~exist( imageFile, 'file')
    fprintf('Cannot load image file %s\n', imageFile);
    return;
else
    oni = imread(imageFile);
    osz = size(oni);
    f1 = nerve_r*2/osz(1);
    onibig = imresize(oni, f1, 'bicubic');
    M.onibigbw = imbinarize(rgb2gray(onibig),'adaptive');
end


%% Init

if ~isempty(file)
    file_full_addr = ['models\' file '.mat'];
    if exist(file_full_addr,'file') && ~any(f('rewrite'))
        fprintf('Reading existing model: %s\n', file);
        load(file_full_addr, 'M');
        if any(f('plot_model')), M.plot.model(); end
        return;
    end
end

% Create bundles
bund_g = [];
bundle_dens = 1;
ok = 0;
ip = 0;

if any(f('GUI'))
    fig_ui = figure('units', 'normalized', 'toolbar','figure', 'Name', 'Bundle Settings', 'outerposition', [0.2 0.1 .6 .8]);
    ax_ui = axes('units', 'normalized', 'position', [0.1 0.1 .7 .8]);
    
    h_feedback = pair_text('Status: ', [.05 .92], []); h_feedback.Position = [.13 .92 .8 .05]; h_feedback.Enable = 'off';
    h_nerv_r = pair_text('Nerve Radius:', [.82 .8], nerve_r);
    h_min_dis = pair_text('Min Clearance', [.82 .7], min_bundles_dis);
    h_bundle_r_range = pair_text('Bundle Raduis Range:', [.82 .6], bundle_r_range);
    %h_neuron_scale = pair_text('Neuron Scale', [.82 .5], neuron_scale);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .4 .12 .05], 'String', 'Gen Bundles', 'Callback', @genBundle);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .3 .12 .05], 'String', 'Gen Neurons', 'Callback', @genNeurons);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .2 .12 .05], 'String', 'Gen Mesh', 'Callback', @genMesh);
    genBundle(0,0);
    
    waitfor(fig_ui);
    if ~ok
        return;
    end
else
    n_tries = 20;
    while isempty(bund_g) && n_tries > 0
        bund_g = fill_circles(0, nerve_r, bundle_dens, bundle_r_range, min_bundles_dis, 7, ~obstaclesOnBundles, []);
        n_tries = n_tries - 1;
        fprintf("Fill Circles: try# %d\n", n_tries);
    end
    if ~n_tries
        warning('Could not fit any bundles within the nerve, please adjust the parameters');
        return;
    end
end

%% Bundle GUI


    function genBundle(~,~)
        t1 = str2num(h_nerv_r.String); %#ok<*ST2NM>
        t2 = str2num(h_min_dis.String);
        t3 = str2num(h_bundle_r_range.String);
        %t4 = str2num(h_neuron_scale.String);
        if isempty(t1) || ~all(size(t1) == [1 1])
            h_feedback.String = 'Bad Nerve Radius!';
            return;
        end
        if isempty(t2) || ~all(size(t2) == [1 1])
            h_feedback.String = 'Bad Bundle Density!';
            return;
        end
        if isempty(t3) || ~all(size(t3) == [1 2])
            h_feedback.String = 'Bad Bundle Radius Range!';
            return;
        end
        nerve_r = t1;
        min_bundles_dis = t2;
        bundle_r_range = t3;
        %neuron_scale = t4;

        n_tries = 20;
        bund_g = [];
        while isempty(bund_g) && n_tries > 0
            bund_g = fill_circles(0, nerve_r, bundle_dens, bundle_r_range, min_bundles_dis, 7, ~obstaclesOnBundles, []);
            n_tries = n_tries - 1;
            fprintf("Fill Circles: try# %d\n", n_tries);
        end
        if ~n_tries
            h_feedback.String = 'Could not fit any bundles within the nerve, please adjust the parameters';
            return;
        end

        osz = size(oni);
        f1 = nerve_r*2/osz(1);
        onibig = imresize(oni, f1, 'bicubic');
        M.onibigbw = imbinarize(rgb2gray(onibig),'adaptive');

        M.bund = bund_g;
        M.nerve_r = nerve_r;
        axis(ax_ui); cla;
        plot_model('no_neuron');
        h_feedback.String = sprintf('Model updated! (%d bundle(s))', size(bund_g, 2));
    end

    function btn_Down(~,~)
        fprintf("button down called\n");
        cp = get(gca,'currentpoint');
        init_insult = cp(1,1:2);
        set(ip, 'XData', init_insult(1) ,'YData', init_insult(2));
    end

    function h = pair_text(txt, pos, def)
        uicontrol('style','edit', 'Units','Normalized', 'position', [pos .08 .05], 'String', txt, 'Enable', 'off', 'BackgroundColor', [0.9400 0.9400 0.9400]);
        h = uicontrol('style','edit', 'Units','Normalized', 'position', [pos+[.08 0] .08 .05], 'String', num2str(def));
    end
     
    function genMesh(~,~)
        h_feedback.String = "Generating Mesh! Please wait ...";
        drawnow;
        err = generate_mesh(refine); % Generate Mesh
        if err 
            h_feedback.String = "Error while generating mesh. Please update the geometry!";
            drawnow;
        else
            h_feedback.String = "Mesh Generation Complete!";
            drawnow;
            plot_mesh();
            hc = allchild(ax_ui);
            for q = 1:length(hc)
                set(hc(q), 'HitTest', 'off');
            end
            ip = plot(init_insult(1), init_insult(2), 'r.', 'markerSize', 20, 'HitTest', 'off');
            set(gca, 'ButtonDownFcn', @btn_Down);
            h_feedback.String = "Please select the starting point of the injury!";
            drawnow;
            uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .1 .12 .05], 'String', 'RunAlgo', 'Callback', @RunAlgo);
        end
    end
    function genNeurons(~,~)
        %%delete(fig_ui);
        ok = 1;
        h_feedback.String = "Generating Neurons! Please wait ...";
        drawnow;
        tic;
        generating_neurons();
        toc;
        plot_model();
        drawnow;
        h_feedback.String = "Generation Done";
        drawnow;

    end

    function RunAlgo(~,~)
        propagation_alg();
        h_feedback.String = "Simulation Animation!";
        M.P.anim('plot_mesh', 'show_dots', 'dt', 200); % playback
        h_feedback.String = "Simulation Complete! Use 'Close' to terminate the window!";
    end

    function runGUI()

    end
%% Mesh Generation
        % finds the average neuron size based on statistics in: Pan et. al [2012]
    function r = radius_avg(norm_cent)
        if ~isfield(M, 'r_dist')
            % Radius average in different regions
            v = [1.04 1.09	1.15    1.21    1.28    1.27
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
        surf(M.r_dist.xq, M.r_dist.yq, M.r_dist.r_avg, 'edgecolor', 'red'); view(0, 90);
        axis equal
    end


    function plot_model(varargin)
        ff = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
        no_neuron = any(ff('no_neuron')); varargin(no_neuron) = [];

        axis equal, hold on
        draw_circles([0 0], M.nerve_r, 200, true, false);
        bundleFill = ~false;
        if ~no_neuron
            draw_circles(M.bund(1:2,:)', M.bund(3,:)', 100, bundleFill, true);
            mlen = length(M.neuron);
            for kk = 1:mlen
                tic;
                draw_circles(M.neuron{kk}(1:2,:)', M.neuron{kk}(3,:)', 5, ~bundleFill, false);
                fprintf("Iter %d out of %d\n", kk, mlen);
                toc;
            end
        else
           draw_circles(M.bund(1:2,:)', M.bund(3,:)', 100, bundleFill, false); 
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
            pos = plotboxpos(par_ax);
            for id = 1:length(M.neuron)
                set(ax(id), 'position', [pos(1) + (p(id,1)-min(xl))/diff(xl)*pos(3),...
                    pos(2) + (p(id,2)-min(yl))/diff(yl)*pos(4),...
                    p(id,3)/diff(xl)*pos(3), p(id,4)/diff(yl)*pos(4)]);
            end
        end
    end


    function generating_neurons

        fprintf('Generating model:\n');

        % whole_geom = [[0 0 nerve_r]' bund_g];
        n_bunds = size(bund_g, 2);
        fprintf('\tNumber of bundles: %d\n', n_bunds);

        neuron_g = cell(n_bunds,1);
        expected_r_avg = zeros(n_bunds,1);

        fprintf('\tNumber of neurons: ');
        total = 0;
        for k = 1:n_bunds
            expected_r_avg(k) = radius_avg(bund_g(1:2,k)./nerve_r);
            %neuron_g{k} = fill_circles(expected_r_avg(k), bund_g(3,k), neuron_dens, neuron_r_range, min(neuron_r_range)/2);
            neuron_g{k} = fill_circles(5, bund_g(3,k), neuron_dens, neuron_r_range, min_bundles_dis/3, 0, ~obstaclesOnBundles, M.onibigbw);
            total = total + size(neuron_g{k}, 2);
            neuron_g{k}(1,:) = neuron_g{k}(1,:) + bund_g(1,k);
            neuron_g{k}(2,:) = neuron_g{k}(2,:) + bund_g(2,k);
        end
        fprintf('Total neurons %d\n', total);

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

    end

    function plot_mesh
        fprintf('Plotting... ');
        h = pdemesh(M.mesh.p,M.mesh.e,M.mesh.t); axis equal
        h(2).Visible = 'off';
        % set(h(2), 'LineWidth', 2, 'Color', [.9 .2 .5]); % Boundaries linewidth and color
        hold on, M.plot.model('LineWidth', 2);
        fprintf('DONE\n');
        drawnow
    end

    function err = generate_mesh(refine)

        err = false;
        fprintf('\tCreating CSG model... ');
        %whole_geom = [[0 0 M.nerve_r]' M.bund M.neuron{:}];
        whole_geom = [[0 0 M.nerve_r]' M.neuron{:}];
        size(whole_geom)
        %whole_geom
        theModel =  [ones(1,size(whole_geom,2)); whole_geom];
        gstat = csgchk(theModel);
        if any(gstat)
           fprintf("Invalid model. Please try again\n"); 
           err = true;
           return;
        end
        try
            
            [M.csg.dl, M.csg.bt] = decsg(theModel);
        catch
            fprintf("Invalid model. Please try again\n"); 
            err = true;
            return;
        end
        fprintf('DONE\n');

        dl = M.csg.dl;

        fprintf('Generating mesh... ');
        try
            [p,e,t] = initmesh(dl);
        catch
            err = true;
            return;
        end
        fprintf('DONE\n');

        if refine > 0
            fprintf('Refining mesh... ');
            for k = 1:refine
                [p,e,t] = refinemesh(dl, p,e,t);
            end
            fprintf('DONE\n');
        end

        p = jigglemesh(p,e,t);

        M.mesh.p = p; M.mesh.e = e; M.mesh.t = t;

        M.plot.mesh = @plot_mesh;
        return;

    end


    function propagation_alg
    % Propagation based on closest points and length of connections

    neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2

    %% Find each domain speed

    bt = M.csg.bt;
    nb_domains = length(bt);
    % Find mapping between domains and geometry elements

    dom_map = [find(sum(bt, 2) == 1);... % Nerve
        find(sum(bt, 2) == 2);... % Bundles
        find(sum(bt, 2) == 3)]; % Neurons

    subdomain.meduim = dom_map(1);
    subdomain.bundles = dom_map(1 + (1:length(M.neuron)));
    all_neurons_g = [M.neuron{:}];
    subdomain.neuron = dom_map(1 + length(M.neuron) + (1:length(all_neurons_g)));

    % assigning domain speed
    domain_speed = zeros(nb_domains,1);
    domain_speed(subdomain.meduim) = medium_speed;
    domain_speed(subdomain.bundles) = bundle_speed;
    domain_speed(subdomain.neuron) = neuron_speed * neuron_speed_formula(all_neurons_g(3,:));

    p = M.mesh.p; e = M.mesh.e; t = M.mesh.t;
    % Setting propagation speed to each point based on their domain speed
    propagation_speed = ones(1, length(p)); % speed of each point
    edge_delay = zeros(1, length(p)); % additional imagenary distance infection has to go in order to pass the boundaries

    for k = 1:nb_domains
        [domain_inter, domain_edge] = pdesdp(p,e,t,k); % find each domain edge and interior points
        propagation_speed(domain_edge) = medium_speed;
        propagation_speed(domain_inter) = domain_speed(k);
    end
    Edge = e(1,:); % set of points on any boundary (slower propagation)
    edge_delay(Edge) = all_edge_delay;

    %% Init

    Conn = t(1:3,:); % Connections / Connectivity matrix
    active_cols_mask = true(length(Conn),1); % columns of Conn matrix that have a dead point
    newly_dead_col_mask = false(length(Conn),1);

    newly_dead = [];
    [~,id] = min(((p(1,:) - init_insult(1)).^2 + (p(2,:) - init_insult(2)).^2)); % id: Closest mesh point to the starting point
    newly_infected = [id; 0];
    infected = []; % first row is the infected point index, second row is the death timer
    death_time = zeros(length(p),1);
    Now = 0; % current time

    %% Propogation algorithm
    % ************************************************************
    % Point states are: live, infected, dead. Those infected have a death timer
    % which will kill them after a certain time. While in infected state A
    % single point can get infected from several sorrounding points, each time,
    % death timer will be updated to the minimum death time assigned.
    % ************************************************************

    fprintf('Simulation... ');

    del = @(x) fprintf(repmat('\b',1,x));
    nDispChar = 0; last_disp_num = 0;

    while any(active_cols_mask) || ~isempty(infected)

        new_disp_num = round(100*(1 - sum(active_cols_mask)/length(Conn)));
        if new_disp_num ~= last_disp_num
            del(nDispChar); nDispChar = fprintf('%.0f%%\t', new_disp_num);
            last_disp_num = new_disp_num;
        end
        if ~isempty(newly_dead)
            % ************************************************************
            % finds all connected points to the newly dead point and based on
            % the length of connection assigns a death timer to those points
            % results go to newly_infected (first row is the infected point
            % index, second row is the death timer)
            % ************************************************************
            newly_dead_col_idx = ceil(find(Conn == newly_dead)/3);

            newly_dead_col_mask(:) = 0;
            newly_dead_col_mask(newly_dead_col_idx) = 1;
            x = Conn(:,active_cols_mask & newly_dead_col_mask); % all surrounding connected points (repetitive and including current point idx)
            y = unique(x(:)); % unique points (including current point idx)
            cpi = y(y ~= newly_dead)'; % connected point index

            dp = sum((p(:,newly_dead)*ones(1,length(cpi)) - p(:,cpi)).^2); % distance of connected points to current dead point
            death_timer = dp./propagation_speed(cpi) + edge_delay(cpi); % death timer is set based on distance, current meduim speed and if that point is on a boundary
            newly_infected = [cpi; death_timer]; % first row is the infected point index, second row is the death timer
            active_cols_mask(newly_dead_col_mask) = 0; % update active connection set

        end
        if ~isempty(infected)
            % ************************************************************
            % Finds points already infected and checks if the new infection
            % death timer is less than the death time already assigned, if so
            % updates the death timer.
            % ************************************************************
            rm = []; % cols to be removed from newly_infected (those already infected)
            for k = 1:size(newly_infected,2)
                idx = find(infected(1,:) == newly_infected(1,k)); % find intersection of two sets
                if ~isempty(idx)
                    infected(2,idx) = min(infected(2,idx), newly_infected(2,k)); % choose minimum death time
                    rm = [rm k];
                end
            end
            newly_infected(:,rm) = []; % remove already-infected points from the set
        end

        infected = [infected newly_infected]; % add non-repetitive newly-infected points to the infected set

        [Dt, idx] = min(infected(2,:)); % find the next infected point waiting for death and its timer value
        Now = Now + Dt; % update time to that tragic time
        newly_dead = infected(1,idx); % mark that point as newly_dead for next loop iteration
        death_time(newly_dead) = Now; % record its absolute death time
        infected(:,idx) = []; % remove point from infected set
        infected(2,:) = infected(2,:) - Dt; % update all death timers
    end
    fprintf('\n');

    end_time = Now;

    P.init_insult = init_insult;
    P.death_time = death_time;
    P.end_time = end_time;
    P.anim = @anim;
    M.P = P;

    %% Play (and record) Animation
        function anim(varargin)
            f = @(f) (cellfun(@(x) ischar(x) && strcmp(x,f), varargin));
            function out = read_pair(name,def)
                ind = find(f(name));
                if ind
                    out = varargin{ind + 1};
                else
                    out = def;
                end
            end
            movie_file = read_pair('movie_file', ''); % Record animation name (records if not empty)
            dt = read_pair('dt', 0.2);
            show_cont = any(f('show_cont')); % show infection contour
            show_dots = any(f('show_dots')); % show infected mesh points
            newfigure = any(f('figure')); % draw in new figure

            record_flag = ~isempty(movie_file); % record flag

            if record_flag
                wO = VideoWriter(movie_file,'MPEG-4');
                wO.FrameRate = 20;
                wO.Quality = 100;
                wO.open;
            end

            if newfigure
                fig = figure('units', 'normalized', 'toolbar','figure', 'Name', 'Bundle Settings', 'outerposition', [0.2 0.15 .6 .8]);
                ax = axes('units', 'normalized', 'position', [0.1 0.1 .8 .8]);
                sld = uicontrol('style','slider', 'Units','Normalized', 'position', [.1 .95 .8 .05]);
            end
            M.plot.model();
            if any(f('plot_mesh'))
                M.plot.mesh();
            end
            plot(init_insult(1), init_insult(2), 'b.', 'markersize', 20);

            if show_cont
                upt = Contour(p, t, death_time);
            end

            if show_dots
                h_new = plot(inf, inf, 'g.', 'markersize', 10);
                h_inf = plot(inf, inf, 'r.', 'markersize', 10);
                setXY = @(h, p) set(h, 'Xdata', p(1,:), 'YData', p(2,:));
            end
            sld.Value = .51;
            tim = 0;
            while true
                try
                    dt = ((sld.Value-.5)*end_time)/20;
                    tim = tim + dt;
                    tim = max(tim, dt);
                    tim = min(tim, end_time);
                    if show_dots
                        setXY(h_new, p(:,tim - dt < death_time & death_time < tim));
                        setXY(h_inf, p(:,tim  > death_time));
                    end
                    if show_cont, upt(tim); end
                    drawnow, if record_flag, wO.writeVideo(getframe); end
                catch
                    break;
                end
            end
            if record_flag, wO.close; end
        end
    end


end


