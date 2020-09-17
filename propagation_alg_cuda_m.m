
function propagation_alg_cuda_m(simIterations, M)
% Propagation based on closest points and length of connections

%neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2
% Load the kernel
filePrefix = 'propagation_algo_cuda';
cudaFilename = strcat(filePrefix, '.cu');
M.h_feedback.String = sprintf('Compiling CUDA file* %s', cudaFilename);
drawnow;
if system(sprintf('nvcc -arch sm_50 -ptx  %s', cudaFilename)) ~= 0
    M.h_feedback.String = 'Error: nvcc compile errors. Aborting Simulation.';
    drawnow;
    return;
end

M.h_feedback.String = sprintf('Compilation successful! Building GPU Data...');
drawnow;

ptxFilename =  strcat(filePrefix, '.ptx');
kernel = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );


kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock,1,1];
kernel.GridSize = [ceil(M.nerve_r*M.nerve_r/kernel.MaxThreadsPerBlock),1];

xinsult = M.init_insult(1);
yinsult = M.init_insult(2);

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


gMap1 = gpuArray(cast(M.cMap1, 'single'));
gMap2 = gpuArray(cast(M.cMap2, 'single'));
gpox  = gpuArray(cast(M.poxMap, 'single'));
gscav = gpuArray(cast(M.scavMap,'single'));
gsmap = gpuArray(cast(M.spaceMap, 'int32'));
%gpixg  = gpuArray(M.pixelMap(:,:,2));
%gpixb  = gpuArray(M.pixelMap(:,:,3));
diff = M.diffInside;

tox_switch = 1;
M.h_feedback.String = 'Running simulation .... ';
drawnow;

for i=1:simIterations/2
    %tic
    gMap2 = feval( kernel, M.nerve_r, tox_switch, gMap1, gMap2, gpox, gscav, gsmap, diff );
    
    tox_switch = ~tox_switch;
    
    gMap1 = feval( kernel, M.nerve_r, tox_switch, gMap1, gMap2, gpox, gscav, gsmap, diff );
    
    %toc
    if mod(i, 50) == 1
        M.h_feedback.String = sprintf('Running simulation: [Iter %d of %d]', i, simIterations);
        %fprintf("==>> %d\n",i);
          rpixel = gather(gMap1);
%         gpixel = gather(gpix);
%         bpixel = gather(bpix);
%         pm = cat(3, rpixel, gpixel, bpixel);
        %                plot_model('pixelmap');
          figure(M.fig_ui);
          axis equal, hold on
          imagesc(rpixel, 'XData', [-1*M.nerve_r M.nerve_r], 'YData', [-1*M.nerve_r, M.nerve_r]);
          axis('on', 'image');
          hold off;
        drawnow
    end
end
M.h_feedback.String = 'Simulation Complete';
end