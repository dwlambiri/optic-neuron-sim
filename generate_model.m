function M = generate_model(varargin)
% GENERATE_MODEL
% parameters:
%   Type, Value
%       neuron_scale
%       nerve_r
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
bundle_r_range = read_pair('bundle_r_range', [nerve_r-1 nerve_r]);
min_bundles_dis = read_pair('min_bundles_dis', 0);
neuron_r_range = [1.5 30]; % because of 10 pixels/um and radius
neuron_dens = read_pair('neuron_dens', 1);
refine = read_pair('refine', 0);
modelFileName = read_pair('file', []);
imageFile = read_pair('image', []);
mielinWidth = read_pair('mielin', 5);
init_insult = read_pair('init_insult', [nerve_r nerve_r]);
insultAmount = read_pair('insult_amount', 25);
insultRadius = read_pair('insult_radius', nerve_r/5);
scavIn = read_pair('scavenging_intra', 0.01);
scavOut = read_pair('scavenging_out', 0.001);
prodAmount = read_pair('production_amount', 0.01);
deathThr   = read_pair('death_threashold', 22);
deathRelease = read_pair('death_release', 10000);
diffusionInside = read_pair('diffusion_inside', 0.02);
diffusionOutside = read_pair('diffusion_outside', 0.02);
diffusionBoundary = read_pair('diffusion_boundary', 0.02);
simIterations = read_pair('iterations', 5000);

if any(f('zoned'))
    zoned_r = true;
else
    zoned_r = false;
end


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

M = struct;

if ~isempty(modelFileName)
    file_full_addr = ['models\' modelFileName];
    if exist(file_full_addr,'file') && ~any(f('rewrite'))
        fprintf('Reading existing model: %s\n', modelFileName);
        load(file_full_addr, 'M');
        if any(f('plot_model')), M.plot.model(); end
        return;
    end
else
    modelFileName = 'temp.mat';
   file_full_addr = ['models\temp.mat'];
end

M.file_full_addr = file_full_addr;
M.save = @(M) save(M.file_full_addr, 'M');

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
    h_mielin = pair_text('Max Mielin Width:', [.82 .55], mielinWidth);
    h_modelFile = pair_text('Model File:', [.82 .5], modelFileName);
    
    %h_neuron_scale = pair_text('Neuron Scale', [.82 .5], neuron_scale);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .4 .12 .05], 'String', 'Load Model', 'Callback', @loadModel);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .35 .12 .05], 'String', 'Save Model', 'Callback', @saveModel);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .3 .12 .05], 'String', 'Gen Bundles', 'Callback', @genBundle);
    uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .25 .12 .05], 'String', 'Gen Neurons', 'Callback', @genNeurons);
    %uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .30 .12 .05], 'String', 'Plot Hist', 'Callback', @runHistogram);
    genBundle(0,0);
    
    waitfor(fig_ui);
    if ~ok
        return;
    end
else
    n_tries = 20;
    while isempty(bund_g) && n_tries > 0
        [bund_g, b] = fill_circles(h_feedback,0, nerve_r, bundle_dens, bundle_r_range, min_bundles_dis, 7, ~obstaclesOnBundles, [], false, 0);
        n_tries = n_tries - 1;
        fprintf("Fill Circles: try# %d\n", n_tries);
    end
    if ~n_tries
        warning('Could not fit any bundles within the nerve, please adjust the parameters');
        return;
    end
end

%% Bundle GUI
    function saveModel(~,~)
        h_feedback = 'Saving Model';
        drawnow;
        t = h_modelFile.String;
        if isempty(t) 
           t = 'tmp.mat'; 
        end
        M.file_full_addr = ['models\' t];
        M.save(M);
        h_feedback = 'Model saved';
        drawnow;
    end

    function loadModel(~,~)
        h_feedback = 'Reading Model';
        drawnow;
        t = h_modelFile.String;
        if isempty(t) 
           t = 'tmp.mat'; 
        end
         M.file_full_addr = ['models\' t];
         if exist(M.file_full_addr,'file')
            h_feedback = sprintf('Reading existing model: %s\n', M.file_full_addr);
            drawnow;
            load(file_full_addr, 'N');
            M.plot.model();
            h_feedback = 'New model loaded'; 
            drawnow;
            return;
         else
            h_feedback = 'Error: Model file not found'; 
            drawnow;
         end
    end

    function genBundle(~,~)
        t1 = str2num(h_nerv_r.String); %#ok<*ST2NM>
        t2 = str2num(h_min_dis.String);
        t3 = str2num(h_bundle_r_range.String);
        t4 = str2num(h_neuron_r_range.String);
        t5 = h_image_file.String;
        t6 = str2num(h_mielin.String);
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
        if isempty(t6)
            t6 = 0;
        end
        nerve_r = t1;
        init_insult = [nerve_r nerve_r];
        min_bundles_dis = t2;
        bundle_r_range = t3;
        neuron_r_range = t4;
        imageFile = t5;
        mielinWidth = t6;
        
        n_tries = 20;
        bund_g = [];
        while isempty(bund_g) && n_tries > 0
            [bund_g, drop] = fill_circles(h_feedback,0, nerve_r, bundle_dens, bundle_r_range, min_bundles_dis, 7, ~obstaclesOnBundles, [],false, 0);
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
        total = generating_neurons();
        %toc;
        plot_model('pixelmap');
        h_feedback.String = sprintf('Generation Done. Model has %d axons. Elapsed time %s seconds.', total, toc);
        drawnow;
        %         drawnow;
        %         figure(M.fig_ui);
        %         ip = plot(init_insult(1), init_insult(2), 'r.', 'markerSize', 20, 'HitTest', 'off');
        %         set(gca, 'ButtonDownFcn', @btn_Down);
        %         h_feedback.String = "Please select the starting point of the injury!";
        %         drawnow;
        uicontrol('style','pushbutton', 'Units','Normalized', 'position', [.84 .1 .12 .05], 'String', 'RunAlgo', 'Callback', @RunAlgo);
        
        
    end

    function RunAlgo(~,~)
        propagation_alg_cuda_m(simIterations, M);
       % h_feedback.String = "Simulation Animation!";
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
        
        figure(M.fig_ui);
        % draw pixel map
        if any(ff('pixelmap'))
            axis equal, hold on
            imagesc(M.pixelMap, 'XData', [-1*nerve_r nerve_r], 'YData', [-1*nerve_r, nerve_r]);
            axis('on', 'image');
            hold off;
            drawnow
            return;
        end
        
        % draw circles
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


    function total = generating_neurons()
        
        h_feedback.String ='Generating model .';
        drawnow;
        
        % whole_geom = [[0 0 nerve_r]' bund_g];
        n_bunds = size(bund_g, 2);
        fprintf('\tNumber of bundles: %d\n', n_bunds);
        
        neuron_g = cell(n_bunds,1);
        blood_g  = cell(n_bunds,1);
        expected_r_avg = zeros(n_bunds,1);
        
        fprintf('\tNumber of neurons: ');
        total = 0;
        for k = 1:n_bunds
            expected_r_avg(k) = radius_avg(bund_g(1:2,k)./nerve_r);
            [neuron_g{k} blood_g{k}] = fill_circles(h_feedback, 5, bund_g(3,k), neuron_dens, neuron_r_range, min_bundles_dis/3, 3, obstaclesOnBundles, M.onibigbw, zoned_r, mielinWidth);
            total = total + size(neuron_g{k}, 2);
            neuron_g{k}(1,:) = neuron_g{k}(1,:) + bund_g(1,k);
            neuron_g{k}(2,:) = neuron_g{k}(2,:) + bund_g(2,k);
            if ~isempty(blood_g{k})
                blood_g{k}(1,:) = blood_g{k}(1,:) + bund_g(1,k);
                blood_g{k}(2,:) = blood_g{k}(2,:) + bund_g(2,k);
            end
        end
        h_feedback.String = sprintf('Model has %d axons', total);
        drawnow;
        
        fprintf('Generating Pixel Map\n');
        pixelMap = cat(3, zeros(2*nerve_r, 2*nerve_r), zeros(2*nerve_r, 2*nerve_r),zeros(2*nerve_r, 2*nerve_r));
        spaceMap = -1*ones(2*nerve_r, 2*nerve_r);
        centerMap = -1*ones(2*nerve_r, 2*nerve_r);
        axonMap = -1*ones(2*nerve_r, 2*nerve_r);
        concentrationMap1 = zeros(2*nerve_r, 2*nerve_r);
        concentrationMap2 = zeros(2*nerve_r, 2*nerve_r);
        poxMap = zeros(2*nerve_r, 2*nerve_r);
        scavMap = zeros(2*nerve_r, 2*nerve_r);
        
        for i = 1: 2*nerve_r
            for j = 1: 2*nerve_r
                if sqrt((i-nerve_r)^2+(j-nerve_r)^2) <= nerve_r
                    %pixelMap(i,j,:) = 0.4;
                    spaceMap(i,j) = 0;
                    centerMap(i, j) = 0;
                end
                if ~isempty(blood_g{k})
                    x = blood_g{k}(1,1);
                    y = blood_g{k}(2,1);
                    r = blood_g{k}(3,1);
                    if sqrt((i-x)^2+(j-y)^2) <= r
                      %pixelMap(i,j,:) = 0.4;
                      spaceMap(i,j) = -1;
                      centerMap(i, j) = -1;
                    end
                    x = blood_g{k}(1,2);
                    y = blood_g{k}(2,2);
                    r = blood_g{k}(3,2);
                    if sqrt((i-x)^2+(j-y)^2) <= r
                      %pixelMap(i,j,:) = 0.4;
                      spaceMap(i,j) = -1;
                      centerMap(i, j) = -1;
                    end
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
                m  = neuron_g{k}(4,q);
                xmin = ceil(xc+nerve_r - (r+m));
                xmax = ceil(xc+nerve_r + (r+m));
                ymin = ceil(nerve_r -yc - (r+m));
                ymax = ceil(nerve_r -yc + (r+m));
                xcn = floor(xc)+nerve_r;
                ycn = nerve_r - floor(yc);
                axonMap(ycn, xcn) = 1;
                for i = ymin:ymax
                    for j = xmin:xmax
                        if sqrt((i-ycn)^2+(j-xcn)^2) <= (r+m)
                            if sqrt((i-ycn)^2+(j-xcn)^2) >= (r-m)
                                pixelMap(i,j,1) = 0;
                                pixelMap(i,j,2) = 0.8;
                                pixelMap(i,j,3) = 0;
                            else
                                pixelMap(i,j,1) = 0;
                                pixelMap(i,j,2) = 1;
                                pixelMap(i,j,3) = 0;
                            end
                            spaceMap(i, j ) = (k-1)*bconst + q; % neuron index
                            linIndex = (xcn-1)*2*nerve_r+ ycn;
                            centerMap(i,j) = linIndex;
%                             if axonMap(linIndex) ~= 1
%                                 fprintf('Error in center %d\n', linIndex);
%                             end
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
                    % the production is proportional to the *perimeter*
                    % which is 2*pi*r. pox*pi*r^2 = 2*prod/r*pi*r^2 =
                    % 2*pi*prod*r ==>> total poxProd per neuron ~2*pi*r
                    poxMap(i,j) = 2*prodAmount/neuron(3);
                    scavMap(i,j) = (1-scavIn);
                else
                    %pixel is part of connecting tissue
                    poxMap(i,j) = 0;
                    scavMap(i,j) = (1-scavOut);
                end
            end
        end
        
        fprintf('Done Pixel Map\n');
        
        M.h_feedback = h_feedback;
        M.nerve_r = nerve_r;
        M.init_insult  = init_insult;
        M.insultAmount = insultAmount;
        M.insultRadius = insultRadius;
        M.neuron_r_range = neuron_r_range;
        M.bund = bund_g;
        M.neuron = neuron_g;
        M.expected_r_avg = expected_r_avg;
        M.bconst = bconst;
        M.pixelMap = pixelMap;
        M.spaceMap = spaceMap;
        M.cMap1 = concentrationMap1;
        M.cMap2 = concentrationMap2;
        M.axonMap = axonMap;
        M.centerMap = centerMap;
        M.poxMap = poxMap;
        M.scavMap = scavMap;
        M.diffInside = diffusionInside;
        M.diffBorder = diffusionBoundary;
        M.diffOutside = diffusionOutside;
        M.deathThr = deathThr;
        M.deathRelease = deathRelease;
        M.scavIn = scavIn;
        M.scavOut = scavOut;
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


    function err = generate_delaunay()
        
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


end


