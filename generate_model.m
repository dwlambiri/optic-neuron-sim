function M = generate_model(varargin)
% GENERATE_MODEL
% parameters:
%   Type, Value
%       opticNerveScale_r
%       opticNerveRadius_r
%       min_bundle_dis
%       axonRadiusRange_r
%       axonDensity_r
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

OPTIC_NERVE_RADIUS_C = 750 % optic nerve radius

opticNerveScale_r = read_pair('optic_nerve_scale', 10);
modelResolution_r = read_pair('resolution', 10); % model resolution in pixels/um
opticNerveRadius_r = ceil(OPTIC_NERVE_RADIUS_C * modelResolution_r * opticNerveScale_r / 100);
bundleRadiusRange_r = read_pair('bundle_radius_range', [opticNerveRadius_r-1 opticNerveRadius_r]);
minBundleDistance_r = read_pair('min_bundle_distance', 0);
axonRadiusRange_r = [0.15 3]; % in um
axonDensity_r = read_pair('axon_density', 1);
refine_r = read_pair('refine', 0);
modelFileName_r = read_pair('model_file', []);
opticalNervePatternImageFile_r = read_pair('pattern_file', []);
mielinWidth_r = read_pair('mielin_width', 5);
init_insult = read_pair('init_insult', [opticNerveRadius_r opticNerveRadius_r]);
insultAmount = read_pair('insult_amount', 0);
insultRadius = read_pair('insult_radius', opticNerveRadius_r/5);
scanInsideAxon_r = read_pair('scav_intra_axon', 0.001);
scavOutsideAxon_r = read_pair('scav_outside_axon', 0.0001);
toxProductionPerArea_r = read_pair('production_amount', 0.0045);
deathToxThreshold_r   = read_pair('death_tox_threashold', 5);
extraToxReleaseOnDeath_r = read_pair('death_tox_release', 0);
diffusionDeadAxon_r = read_pair('diffusion_dead', 0.02);
diffusionInsideAxon_r = read_pair('diffusion_inside', 0.02);
diffusionOutsideAxon_r = read_pair('diffusion_outside', 0.02);
diffusionAxonBoundary_r = read_pair('diffusion_boundary', 0.02);
simIterations_r = read_pair('iterations', 5000);
deathToxVar = read_pair('death_tox_var', 0.4);

if any(f('zoned'))
    zoned_r = true;
else
    zoned_r = false;
end


obstaclesOnBundles = true;

oni = [];


if isempty(opticalNervePatternImageFile_r)
    opticalNervePatternImageFile_r = 'optic-nerve.png';
end

if ~exist( opticalNervePatternImageFile_r, 'file')
    fprintf('Cannot load image file %s\n', opticalNervePatternImageFile_r);
    M.onibigbw = [];
    %return;
else
    oni = imread(opticalNervePatternImageFile_r);
    osz = size(oni);
    f1 = opticNerveRadius_r*2/osz(1);
    onibig = imresize(oni, f1, 'bicubic');
    M.onibigbw = imbinarize(rgb2gray(onibig),'adaptive');
end


%% Init

M = struct;

if ~isempty(modelFileName_r)
    file_full_addr = ['models\' modelFileName_r];
    if exist(file_full_addr,'file') && ~any(f('rewrite'))
        fprintf('Reading existing model: %s\n', modelFileName_r);
        load(file_full_addr, 'M');
        if any(f('plot_model')), M.plot.model(); end
        return;
    end
else
    modelFileName_r = 'temp.mat';
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
    
    
    h_nerv_r = pair_text('Radius 750 um :', [.82 .9], []); h_nerv_r.Enable = 'off';
    h_nerv_r = sprintf('%d pix', opticNerveRadius_r);
    
    h_opticNerveScale_r = pair_text('Scale (%)', [.82 .85], opticNerveScale_r);
    h_model_res  = pair_text('Res (pix/um)', [.82 .8], modelResolution_r);
    
    h_min_dis = pair_text('Min Clearance', [.82 .75], minBundleDistance_r);
    h_bundleRadiusRange_r = pair_text('Bundle Radius Range:', [.82 .7], bundleRadiusRange_r);
    h_axonRadiusRange_r = pair_text('Axon Range (um):', [.82 .65], axonRadiusRange_r);
    h_image_file = pair_text('Pattern File:', [.82 .6], opticalNervePatternImageFile_r);
    h_mielin = pair_text('Max Mielin Width:', [.82 .55], mielinWidth_r);
    h_modelFile = pair_text('Model File:', [.82 .5], modelFileName_r);
    
    
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
        [bund_g, b] = fill_circles(h_feedback,0, opticNerveRadius_r, bundle_dens, bundleRadiusRange_r, minBundleDistance_r, 7, ~obstaclesOnBundles, [], false, 0);
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
        
        t2 = str2num(h_min_dis.String);
        t4 = str2num(h_axonRadiusRange_r.String);
        t5 = h_image_file.String;
        t6 = str2num(h_mielin.String);
        t7 = str2num(h_opticNerveScale_r.String);
        t8 = str2num(h_model_res.String);
        
        if isempty(t2) || ~all(size(t2) == [1 1])
            h_feedback.String = 'Bad Bundle Density!';
            return;
        end
        if isempty(t6)
            t6 = 0;
        end
        
        opticNerveRadius_r = ceil(OPTIC_NERVE_RADIUS_C * modelResolution_r * opticNerveScale_r / 100);
        
        h_nerv_r = sprintf('%d pix', opticNerveRadius_r);
        
        init_insult = [opticNerveRadius_r opticNerveRadius_r];
        minBundleDistance_r = t2;
        bundleRadiusRange_r = [opticNerveRadius_r-1 opticNerveRadius_r];
        h_bundleRadiusRange_r = sprintf('[%d %d]', opticNerveRadius_r-1, opticNerveRadius_r);
        axonRadiusRange_r = t4;
        opticalNervePatternImageFile_r = t5;
        mielinWidth_r = t6;
        opticNerveScale_r = t7;
        modelResolution_r = t8;
        n_tries = 20;
        bund_g = [];
        
        drawnow;

        while isempty(bund_g) && n_tries > 0
            [bund_g, drop] = fill_circles(h_feedback,0, opticNerveRadius_r, bundle_dens, bundleRadiusRange_r, minBundleDistance_r, 7, ~obstaclesOnBundles, [],false, 0);
            n_tries = n_tries - 1;
            fprintf("Fill Circles: try# %d\n", n_tries);
        end
        if ~n_tries
            h_feedback.String = 'Could not fit any bundles within the nerve, please adjust the parameters';
            return;
        end
        
        fprintf("IMAGE FILE = %s\n", t5);
        if ~exist( opticalNervePatternImageFile_r, 'file')
            fprintf('Cannot load image file %s\n', opticalNervePatternImageFile_r);
            M.onibigbw = [];
            %return;
        else
            oni = imread(opticalNervePatternImageFile_r);
            osz = size(oni);
            f1 = opticNerveRadius_r*2/osz(1);
            onibig = imresize(oni, f1, 'bicubic');
            M.onibigbw = imbinarize(rgb2gray(onibig),'adaptive');
        end
        M.bund = bund_g;
        M.opticNerveRadius_r = opticNerveRadius_r;
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
        propagation_alg_cuda_m(simIterations_r, M);
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
                1.01 1.1    1.1     1.22    1.22    1.21  ] * opticNerveScale_r;
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
            imagesc(M.pixelMap, 'XData', [-1*opticNerveRadius_r opticNerveRadius_r], 'YData', [-1*opticNerveRadius_r, opticNerveRadius_r]);
            axis('on', 'image');
            hold off;
            drawnow
            return;
        end
        
        % draw circles
        axis equal, hold on
        draw_circles(M.fig_ui, [0 0], M.opticNerveRadius_r, 200, true, false);
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
                xlim(M.axonRadiusRange_r);
                % set(ax(idx), 'XTickLabel', '');
                set(ax(idx), 'YTickLabel', '');
                text(.3, 1.1, sprintf('%.1f(%.1f)', mean(M.neuron{idx}(3,:)), M.expected_r_avg(idx)), 'Units', 'Normalized', 'Color', [.8 .2 .3]);
            end
            set(fig, 'SizeChangedFcn', @size_changed);
            size_changed(0, 0);
        else
            histogram(M.neuron{varargin{1}}(3,:), opts{:});
            xlim(M.axonRadiusRange_r);
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
        
        % whole_geom = [[0 0 opticNerveRadius_r]' bund_g];
        n_bunds = size(bund_g, 2);
        fprintf('\tNumber of bundles: %d\n', n_bunds);
        
        neuron_g = cell(n_bunds,1);
        blood_g  = cell(n_bunds,1);
        expected_r_avg = zeros(n_bunds,1);
        
        axonPixRange = axonRadiusRange_r * modelResolution_r;
        
        fprintf('\tNumber of neurons: ');
        total = 0;
        for k = 1:n_bunds
            expected_r_avg(k) = radius_avg(bund_g(1:2,k)./opticNerveRadius_r);
            [neuron_g{k} blood_g{k}] = fill_circles(h_feedback, 0.5*modelResolution_r, bund_g(3,k), axonDensity_r, axonPixRange, minBundleDistance_r/3, 3, obstaclesOnBundles, M.onibigbw, zoned_r, mielinWidth_r, OPTIC_NERVE_RADIUS_C*modelResolution_r);
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
        pixelMap = cat(3, zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r), zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r),zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r));
        spaceMap = -1*ones(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        centerMap = -1*ones(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        diffMap = zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        axonMap = -1*ones(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        concentrationMap1 = zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        axonDeathValue = zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        poxMap = zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        scavMap = zeros(2*opticNerveRadius_r, 2*opticNerveRadius_r);
        
        diffValues = [ -1 diffusionDeadAxon_r *(modelResolution_r/10)^2 diffusionInsideAxon_r *(modelResolution_r/10)^2 diffusionAxonBoundary_r *(modelResolution_r/10)^2 diffusionOutsideAxon_r *(modelResolution_r/10)^2];
        diffDead = 2;
        diffInside = 3;
        diffBorder = 4;
        diffOutside = 5;
        
        numBits = 3;
        
        for i = 1: 2*opticNerveRadius_r
            for j = 1: 2*opticNerveRadius_r
                if sqrt((i-opticNerveRadius_r)^2+(j-opticNerveRadius_r)^2) <= opticNerveRadius_r
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
        %temp = opticNerveRadius_r / sqrt(2);
        %lenAdjust = 1;
        for k = 1:n_bunds
            lb  = size(neuron_g{k}, 2);
            for q =1: lb
                xc = neuron_g{k}(1,q);
                yc = neuron_g{k}(2,q);
                r  = neuron_g{k}(3,q);
                m  = neuron_g{k}(4,q);
                xmin = ceil(xc+opticNerveRadius_r - (r+m));
                xmax = ceil(xc+opticNerveRadius_r + (r+m));
                ymin = ceil(opticNerveRadius_r -yc - (r+m));
                ymax = ceil(opticNerveRadius_r -yc + (r+m));
                xcn = floor(xc)+opticNerveRadius_r;
                ycn = opticNerveRadius_r - floor(yc);
                axonMap(ycn, xcn) = 1;
                axonDeathValue(ycn, xcn) = deathToxThreshold_r*(1+deathToxVar*xcn/(2*opticNerveRadius_r))/(modelResolution_r^2);
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
                            linIndex = (xcn-1)*2*opticNerveRadius_r+ ycn;
                            centerMap(i,j) = linIndex;
%                             if axonMap(linIndex) ~= 1
%                                 fprintf('Error in center %d\n', linIndex);
%                             end
                        end
                    end
                end
            end
        end
        
        for i = 1: 2*opticNerveRadius_r
            for j = 1: 2*opticNerveRadius_r
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
                    poxMap(i,j) = 2*toxProductionPerArea_r/(neuron(3)/modelResolution_r)/(modelResolution_r^2);
                    scavMap(i,j) = (1-scanInsideAxon_r);
                else
                    %pixel is part of connecting tissue
                    poxMap(i,j) = 0;
                    scavMap(i,j) = (1-scavOutsideAxon_r);
                end
            end
        end
        
        fprintf('Done Pixel Map\n');
        
        M.h_feedback = h_feedback;
        M.opticNerveReal = OPTIC_NERVE_RADIUS_C;
        M.opticNerveRadiusPixels = opticNerveRadius_r;
        M.opticNerveScale = opticNerveScale_r;
        M.modelResolution = modelResolution_r;
        M.init_insult  = init_insult;
        M.insultAmount = insultAmount;
        M.insultRadius = insultRadius;
        M.axonRadiusRange_r = axonRadiusRange_r;
        M.bund = bund_g;
        M.neuron = neuron_g;
        M.expected_r_avg = expected_r_avg;
        M.bconst = bconst;
        M.pixelMap = pixelMap;
        M.spaceMap = spaceMap;
        M.cMap1 = concentrationMap1;
        M.axonMap = axonMap;
        M.centerMap = centerMap;
        M.poxMap = poxMap;
        M.scavMap = scavMap;
        M.axonDeathValue = axonDeathValue;
        M.diffDead = diffDead;
        M.diffInside = diffInside;
        M.diffBorder = diffBorder;
        M.diffOutside = diffOutside;
        M.diffMap = diffMap;
        M.diffValues = diffValues;
        M.deathToxThreshold_r = deathToxThreshold_r / modelResolution_r^2;
        M.extraToxReleaseOnDeath_r = extraToxReleaseOnDeath_r / modelResolution_r^2;
        M.scanInsideAxon_r = scanInsideAxon_r;
        M.scavOutsideAxon_r = scavOutsideAxon_r;
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
        
        %whole_geom = [[0 0 M.opticNerveRadius_r]' M.bund M.neuron{:}];      
        %ang = (linspace(0,2*pi,100)); % min(round(.5*radius), 15)
        %xp = opticNerveRadius_r*cos(ang);
        %yp = opticNerveRadius_r*sin(ang);
        
        whole_x =  M.neuron{:}(1,:)';
        whole_y =  M.neuron{:}(2,:)';
        tri = delaunay(whole_x,whole_y);
        
        M.whole_x = whole_x;
        M.whole_y = whole_y;
        M.tri = tri;
    end


end


