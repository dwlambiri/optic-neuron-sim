
function propagation_alg_arrayfun(simIterations, M)
% Propagation based on closest points and length of connections

%neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2

xinsult = M.init_insult(1);
yinsult = M.init_insult(2);

M.h_feedback.String = 'Running simulation';
drawnow;

for i=1:2*M.insultRadius
    for j=1:2*M.insultRadius
        if sqrt((i-M.insultRadius)^2+(j-M.insultRadius)^2) > M.insultRadius
            continue;
        end
        M.cMap1(yinsult-M.insultRadius+i, xinsult-M.insultRadius+j) = M.insultAmount;
        M.pixelMap(yinsult-M.insultRadius+i, xinsult-M.insultRadius+j, 1) = 1;
    end
end

if M.spaceMap(yinsult,xinsult) < 0
    M.h_feedback.String = sprintf('Error: Insult position is outside the optical nerve x=%d, y=%d\n', xinsult, yinsult);
    drawnow;
    return;
end


gMap1 = gpuArray(M.cMap1);
gu    = gpuArray(circshift(M.cMap1,-1,1));
gd    = gpuArray(circshift(M.cMap1,1,1));
gl    = gpuArray(circshift(M.cMap1,-1,2));
gr    = gpuArray(circshift(M.cMap1,1,2));
gpox  = gpuArray(M.poxMap);
gscav = gpuArray(M.scavMap);
gsmap = gpuArray(M.spaceMap);
%gpixg  = gpuArray(M.pixelMap(:,:,2));
%gpixb  = gpuArray(M.pixelMap(:,:,3));
diff = ones(2*M.nerve_r, 2*M.nerve_r)* M.diffInside;
di = gpuArray(diff);

for i=1:simIterations
    %tic
    [gMap1, gpox, rpix, gpix, bpix] = arrayfun(@gpu_diffusion, ...
        gMap1, gu, gd, gl, gr, gpox, gscav, gsmap, di);
    gu    = gpuArray(circshift(gMap1,-1,1));
    gd    = gpuArray(circshift(gMap1,1,1));
    gl    = gpuArray(circshift(gMap1,-1,2));
    gr    = gpuArray(circshift(gMap1,1,2));
    
    %toc
    if mod(i, 200) == 1
        M.h_feedback.String = sprintf('Running simulation [Iter %d of %d]', i, simIterations);
        %fprintf("==>> %d\n",i);
        rpixel = gather(gMap1);
        gpixel = gather(gpix);
        bpixel = gather(bpix);
        pm = cat(3, rpixel, gpixel, bpixel);
        %                plot_model('pixelmap');
        figure(M.fig_ui);
        axis equal, hold on
        imagesc(pm, 'XData', [-1*M.nerve_r M.nerve_r], 'YData', [-1*M.nerve_r, M.nerve_r]);
        axis('on', 'image');
        hold off;
        drawnow
    end
end
M.h_feedback.String = 'Simulation Complete';
end