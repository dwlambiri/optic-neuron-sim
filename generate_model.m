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
neuron_r_range = [1.5 35]; % because of 10 pixels/um and radius
neuron_dens = read_pair('neuron_dens', 1);
refine = read_pair('refine', 0);
file = read_pair('file', []);
imageFile = read_pair('image', []);

if any(f('zoned'))
    zoned_r = true;
else
    zoned_r = false;
end

% PROPAGATION ALGO PARAMETERS
init_insult = read_pair('init_insult', [700 700]);
medium_speed = read_pair('medium_speed', .02);
bundle_speed = read_pair('bundle_speed', .02);
neuron_speed = read_pair('neuron_speed', 1);
all_edge_delay = read_pair('edge_delay', 2);

obstaclesOnBundles = true;

oni = []; 


if isempty(imageFile) 
     imageFile = 'optic-nerve.png';
end

if ~exist( imageFile, 'file')
    fprintf('Cannot load image file %s\n', imageFile);
    M.onibigbw = [];
    %return;
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
    if zoned_r
        title_r = 'Bundle Settings (ZONED)';
    else
        title_r = 'Bundle Settings (SINGLE)';
    end
    fig_ui = figure('units', 'normalized', 'toolbar','figure', 'Name', title_r, 'outerposition', [0.2 0.1 .6 .8]);
    ax_ui = axes('units', 'normalized', 'position', [0.1 0.1 .7 .8]);
    
    M.fig_ui = fig_ui;
    M.ax_ui = ax_ui;
    
    h_feedback = pair_text('Status: ', [.05 .92], []); h_feedback.Position = [.13 .92 .8 .05]; h_feedback.Enable = 'off';
    h_nerv_r = pair_text('Nerve Radius:', [.82 .8], nerve_r);
    h_min_dis = pair_text('Min Clearance', [.82 .75], min_bundles_dis);
    h_bundle_r_range = pair_text('Bundle Radius Range:', [.82 .7], bundle_r_range);
    h_neuron_r_range = pair_text('Neuron Radius:', [.82 .65], neuron_r_range);
    h_image_file = pair_text('Image File:', [.82 .6], imageFile);
    %h_neuron_scale = pair_text('Neuron Scale', [.82 .5], neuron_scale);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .4 .12 .05], 'String', 'Gen Bundles', 'Callback', @genBundle);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .35 .12 .05], 'String', 'Gen Neurons', 'Callback', @genNeurons);
    %uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .30 .12 .05], 'String', 'Plot Hist', 'Callback', @runHistogram);
    genBundle(0,0);
    
    waitfor(fig_ui);
    if ~ok
        return;
    end
else
    n_tries = 20;
    while isempty(bund_g) && n_tries > 0
        bund_g = fill_circles(0, nerve_r, bundle_dens, bundle_r_range, min_bundles_dis, 7, ~obstaclesOnBundles, [], false);
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
        t4 = str2num(h_neuron_r_range.String);
        t5 = h_image_file.String;
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
        neuron_r_range = t4;
        imageFile = t5;

        n_tries = 20;
        bund_g = [];
        while isempty(bund_g) && n_tries > 0
            bund_g = fill_circles(0, nerve_r, bundle_dens, bundle_r_range, min_bundles_dis, 7, ~obstaclesOnBundles, [],false);
            n_tries = n_tries - 1;
            fprintf("Fill Circles: try# %d\n", n_tries);
        end
        if ~n_tries
            h_feedback.String = 'Could not fit any bundles within the nerve, please adjust the parameters';
            return;
        end

        fprintf("IMAGE FILE = %s\n", t5);
        if ~exist( imageFile, 'file')
            fprintf('Cannot load image file %s\n', imageFile);
            M.onibigbw = [];
            %return;
        else
            oni = imread(imageFile);
            osz = size(oni);
            f1 = nerve_r*2/osz(1);
            onibig = imresize(oni, f1, 'bicubic');
            M.onibigbw = imbinarize(rgb2gray(onibig),'adaptive');
        end
        M.bund = bund_g;
        M.nerve_r = nerve_r;
        axis(M.ax_ui); cla;
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
        err = generate_delaunay(refine); % Generate Mesh
        if err 
            h_feedback.String = "Error while generating mesh. Please update the geometry!";
            drawnow;
        else
            h_feedback.String = "Mesh Generation Complete!";
            drawnow;
            plot_mesh();
            hc = allchild(M.ax_ui);
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
        %plot_model();
        figure(M.fig_ui);
        imshow(M.spaceMap, 'XData', [-1500 1500], 'YData', [-1500, 1500]);
        drawnow;
        h_feedback.String = "Generation Done";
%         drawnow;
%         figure(M.fig_ui);
%         ip = plot(init_insult(1), init_insult(2), 'r.', 'markerSize', 20, 'HitTest', 'off');
%         set(gca, 'ButtonDownFcn', @btn_Down);
%         h_feedback.String = "Please select the starting point of the injury!";
%         drawnow;
        uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .1 .12 .05], 'String', 'RunAlgo', 'Callback', @RunAlgo);


    end

    function RunAlgo(~,~)
        propagation_alg();
        h_feedback.String = "Simulation Animation!";
    end

    function runHistogram(~,~)
        %propagation_alg();
        h_feedback.String = "Plot Radius Histogram";
        plot_histogram(1);
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
        draw_circles(M.fig_ui, [0 0], M.nerve_r, 200, true, false);
        bundleFill = ~false;
        if ~no_neuron
            draw_circles(M.fig_ui,M.bund(1:2,:)', M.bund(3,:)', 100, bundleFill, true);
            mlen = length(M.neuron);
            for kk = 1:mlen
                tic;
                draw_circles(M.fig_ui, M.neuron{kk}(1:2,:)', M.neuron{kk}(3,:)', 5, ~bundleFill, false);
                fprintf("Iter %d out of %d\n", kk, mlen);
                toc;
            end
        else
           draw_circles(M.fig_ui, M.bund(1:2,:)', M.bund(3,:)', 100, bundleFill, false); 
        end
        hold off;
        drawnow
        
    end

    function plot_histogram(varargin)
        opts = {'Normalization', 'probability'};
        if isempty(varargin)
            fig = figure;
            title('Histogram: actual mean (expected mean)');
            %plot_model('no_neuron');
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
            neuron_g{k} = fill_circles(5, bund_g(3,k), neuron_dens, neuron_r_range, min_bundles_dis/3, 3, obstaclesOnBundles, M.onibigbw, zoned_r);
            total = total + size(neuron_g{k}, 2);
            neuron_g{k}(1,:) = neuron_g{k}(1,:) + bund_g(1,k);
            neuron_g{k}(2,:) = neuron_g{k}(2,:) + bund_g(2,k);
        end
        fprintf('Total neurons %d\n', total);
        
        fprintf('Generating Pixel Map\n');
        pixelMap = 32*ones(2*nerve_r, 2*nerve_r);
        spaceMap = -1*ones(2*nerve_r, 2*nerve_r);
        concentrationMap1 = zeros(2*nerve_r, 2*nerve_r);
        concentrationMap2 = zeros(2*nerve_r, 2*nerve_r);
        poxMap = zeros(2*nerve_r, 2*nerve_r);
        scavMap = zeros(2*nerve_r, 2*nerve_r);
        
        for i = 1: 2*nerve_r
            for j = 1: 2*nerve_r
                if sqrt((i-nerve_r)^2+(j-nerve_r)^2) <= nerve_r
                    pixelMap(i,j) = 0;
                    spaceMap(i,j) = 0;
                end
            end
        end

        bconst = 500000; % [DWL] need some formula based on bundles; large enough for now.
        for k = 1:n_bunds
            lb  = size(neuron_g{k}, 2);
            for q =1: lb
                xc = neuron_g{k}(1,q);
                yc = neuron_g{k}(2,q);
                r  = neuron_g{k}(3,q);
                xmin = ceil(xc+nerve_r - r);
                xmax = ceil(xc+nerve_r + r);
                ymin = ceil(nerve_r -yc - r);
                ymax = ceil(nerve_r -yc + r);
                xcn = ceil(xc+nerve_r);
                ycn = ceil(nerve_r - yc);
                for i = ymin:ymax
                    for j = xmin:xmax
                        if sqrt((i-ycn)^2+(j-xcn)^2) <= r
                            pixelMap(i,j) = 255;
                            spaceMap(i, j ) = (k-1)*bconst + q; % neuron index
                        end
                    end
                end              
            end
        end
        
        for i = 1: 2*nerve_r
            for j = 1: 2*nerve_r
                if spaceMap(i,j) < 0 
                    continue
                end
                if spaceMap(i,j) > 0
                    % pixel is part of neuron
                    qqq = floor(spaceMap(i,j)/bconst)+1;
                    idx = mod(spaceMap(i,j), bconst);
                    neuron = neuron_g{qqq}(:,idx);
                    
                    if neuron(3) < 15
                        poxMap(i,j) = 1;
                        scavMap(i,j) = -0.2;
                    else
                        poxMap(i,j) = 1;
                        scavMap(i,j) = -2;
                    end
                else
                    %pixel is part of connecting tissue
                    poxMap(i,j) = 0;
                    scavMap(i,j) = -10^-2;
                end
            end
        end
        
        fprintf('Done Pixel Map\n');
        
        M.nerve_r = nerve_r;
        M.neuron_r_range = neuron_r_range;
        M.bund = bund_g;
        M.neuron = neuron_g;
        M.expected_r_avg = expected_r_avg;
        M.bconst = bconst;
        M.pixelMap = pixelMap;
        M.spaceMap = spaceMap;
        M.cMap1 = concentrationMap1;
        M.cMap2 = concentrationMap2;
        M.poxMap = poxMap;
        M.scavMap = scavMap;
        M.diffInside = 0.2;
        M.diffBorder = 0.02;
        M.diffOutside = 0.2;
        M.create_csg = @create_csg;
        M.plot.model = @plot_model;
        M.plot.histogram = @plot_histogram;
        M.plot.avg_r = @plot_avg_r;

        M.file_full_addr = file_full_addr;
        M.save = @(M) save(M.file_full_addr, 'M');

    end

    function plot_mesh
        fprintf('Plotting... ');
        triplot(M.tri,M.whole_x,M.whole_y);
        %hold on, M.plot.model('LineWidth', 2);
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


    function err = generate_delaunay(refine)

        err = false;
        fprintf('\tCreating Delaunay mesh... ');
        %whole_geom = [[0 0 M.nerve_r]' M.bund M.neuron{:}];
        
        
        %ang = (linspace(0,2*pi,100)); % min(round(.5*radius), 15)
        %xp = nerve_r*cos(ang);
        %yp = nerve_r*sin(ang);
        whole_x =  M.neuron{:}(1,:)';
        whole_y =  M.neuron{:}(2,:)';
        tri = delaunay(whole_x,whole_y);
        
        M.whole_x = whole_x;
        M.whole_y = whole_y;
        M.tri = tri;
    end


    function propagation_alg()
        % Propagation based on closest points and length of connections
        
        %neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2
        
        init_insult(1)
        init_insult(2)
        
        xinsult = 1500;
        yinsult = 1500;
        
        M.cMap1(yinsult, xinsult) = 100;
        M.pixelMap = zeros(2*M.nerve_r, 2*M.nerve_r);
        M.pixelMap(yinsult, xinsult) = 1;
        
        cmin = 4;
        
        xmin = max(xinsult, cmin);
        xmax = min(xinsult, 2*M.nerve_r - cmin);
        ymin = max(yinsult, cmin);
        ymax = min(yinsult, 2*M.nerve_r - cmin);
        
        if M.spaceMap(yinsult,xinsult) < 0
           fprintf('Insult position is outside the optical nerve x=%d, y=%d\n', xinsult, yinsult);
           return;
        end
        
        for frame = 1: 1480
            xminconst = false;
            xmaxconst = false;
            anyminy = false;
            anymaxy = false;
            for x = xmin-1: xmax+1
                linechange = false;
                for y = ymin-1: ymax+1
%                     if M.spaceMap(y,x) < 0
%                         continue;
%                     end
                    diffusion = M.diffInside*(M.cMap1(y-1,x)-M.cMap1(y,x)) +  M.diffInside*(M.cMap1(y+1,x)-M.cMap1(y,x)) +  M.diffInside*(M.cMap1(y,x-1)-M.cMap1(y,x)) + M.diffInside*(M.cMap1(y,x+1)-M.cMap1(y,x));
                    if M.cMap1(y,x) == 0 && diffusion == 0
                        continue;
                    end
                    linechange = true;
                    if y == (ymin -1) 
                        anyminy = true;
                    end
                    if y == (ymax+1)
                        anymaxy = true;
                    end
                    M.cMap2(y,x) = max(0,M.cMap1(y,x)+ M.poxMap(y,x) + M.scavMap(y,x) + diffusion) ;
                end
                if x == xmin-1 && linechange == false
                    xminconst = true;
                end
                if x == xmax-1 && linechange == false
                    xmaxconst = true;
                end
            end
            if xminconst == true
                xmin = xmin +1;
            end
            if xmaxconst == true
                xmax = xmax - 1;
            end
            
            if anyminy == false
                ymin = ymin +1;
            end
            
            if anymaxy == false
                ymax = ymax -1;
            end
            
            xmin = max(xmin -1, cmin);
            xmax = min(xmax+1, 2*M.nerve_r - cmin);
            ymin = max(ymin-1,  cmin);
            ymax = min(ymax+1, 2*M.nerve_r- cmin);
            
            
            fprintf('frame %d [1] xmin=%d, xmax=%d, ymin=%d, ymax= %d\n', frame, xmin, xmax, ymin, ymax);
            
            xminconst = false;
            xmaxconst = false;
            anyminy = false;
            anymaxy = false;
            for x = xmin-1: xmax+1
                for y = ymin-1: ymax+1
%                     if M.spaceMap(y,x) < 0
%                         continue;
%                     end
                    diffusion = M.diffInside*(M.cMap2(y-1,x)-M.cMap2(y,x)) +  M.diffInside*(M.cMap2(y+1,x)-M.cMap2(y,x)) +  M.diffInside*(M.cMap2(y,x-1)-M.cMap2(y,x)) + M.diffInside*(M.cMap2(y,x+1)-M.cMap2(y,x));
                    if M.cMap2(y,x) == 0 && diffusion == 0
                        continue;
                    end
                    if y == (ymin -1) 
                        anyminy = true;
                    end
                    if y == (ymax+1)
                        anymaxy = true;
                    end
                    M.cMap1(y,x) = max(0,M.cMap2(y,x)+ M.poxMap(y,x) + M.scavMap(y,x) + diffusion) ;
                end
                if x == xmin-1 && linechange == false
                    xminconst = true;
                end
                if x == xmax-1 && linechange == false
                    xmaxconst = true;
                end
            end
            
            if xminconst == true
                xmin = xmin +1;
            end
            if xmaxconst == true
                xmax = xmax - 1;
            end
            if anyminy == false
                ymin = ymin +1;
            end
            
            if anymaxy == false
                ymax = ymax -1;
            end
            xmin = max(xmin -1, cmin);
            xmax = min(xmax+1, 2*M.nerve_r - cmin);
            ymin = max(ymin-1, cmin);
            ymax = min(ymax+1, 2*M.nerve_r- cmin);
            
            fprintf('frame %d [2] xmin=%d, xmax=%d, ymin=%d, ymax= %d\n', frame, xmin, xmax, ymin, ymax);
            
            if mod(frame, 40) == 1
                figure(M.fig_ui);
                image(M.cMap1, 'XData', [-1500 1500], 'YData', [-1500, 1500]);                 
                colorbar('East');
                drawnow;
            end
            
%             for x = 2: 2*M.nerve_r-1
%                 for y = 2: 2*M.nerve_r-1
%                     if M.spaceMap(y,x) < 0
%                         continue;
%                     end
%                     if M.cMap1(y,x) == 0
%                         continue;
%                     end
%                     M.cMap2(y,x) = M.cMap1(y,x)+ M.poxMap(y,x) + M.scavMap(y,x) + M.diffInside*(M.cMap1(y-1,x)-M.cMap1(y,x)) +  M.diffInside*(M.cMap1(y+1,x)-M.cMap1(y,x)) +  M.diffInside*(M.cMap1(y,x-1)-M.cMap1(y,x)) + M.diffInside*(M.cMap1(y,x+1)-M.cMap1(y,x)) ;
%                 end
%             end

%             for x = 2: 2*M.nerve_r-1
%                 for y = 2: 2*M.nerve_r-1
%                     if M.spaceMap(y,x) < 0
%                         continue;
%                     end
%                     if M.cMap1(y,x) == 0
%                         continue;
%                     end
%                     M.cMap1(y,x) = M.cMap2(y,x)+ M.poxMap(y,x) + M.scavMap(y,x) + M.diffInside*(M.cMap2(y-1,x)-M.cMap2(y,x)) +  M.diffInside*(M.cMap2(y+1,x)-M.cMap2(y,x)) +  M.diffInside*(M.cMap2(y,x-1)-M.cMap2(y,x)) + M.diffInside*(M.cMap2(y,x+1)-M.cMap2(y,x)) ;
%                 end
%             end

        end
    end

end


