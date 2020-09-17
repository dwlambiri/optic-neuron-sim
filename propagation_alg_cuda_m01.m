
function propagation_alg_cuda_m01(simIterations, M)

gd = gpuDevice();
reset(gd);
wait(gd);

% Propagation based on closest points and length of connections

%neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2
% Load the kernel
filePrefix = 'propagation_algo_cuda';
cudaFilename = strcat(filePrefix, '.cu');
fprintf('Compiling CUDA file* %s\n', cudaFilename);

if system(sprintf('nvcc -arch sm_50 -ptx  %s', cudaFilename)) ~= 0
    fprintf('Error: nvcc compile errors. Aborting Simulation.\n');
    return;
end

fprintf('Compilation successful! Building GPU Data...\n');


ptxFilename =  strcat(filePrefix, '.ptx');
kernel = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );

threadPerBlock = 32;
N = 2*M.nerve_r;
kernel.ThreadBlockSize = [ threadPerBlock, threadPerBlock,1];
kernel.GridSize = [ceil(N/threadPerBlock), ceil(N/threadPerBlock),1 ];

xinsult = M.init_insult(1);
yinsult = M.init_insult(2);

for i=1:2*M.insultRadius
    for j=1:2*M.insultRadius
        if sqrt((i-M.insultRadius)^2+(j-M.insultRadius)^2) > M.insultRadius
            continue;
        end
        M.cMap1(yinsult-M.insultRadius+i, xinsult-M.insultRadius+j) = 400*M.insultAmount;
        M.pixelMap(yinsult-M.insultRadius+i, xinsult-M.insultRadius+j, 1) = 1;
    end
end

if M.spaceMap(yinsult,xinsult) < 0
    fprintf('Error: Insult position is outside the optical nerve x=%d, y=%d\n', xinsult, yinsult);
    return;
end


gMap1 = gpuArray(cast(M.cMap1, 'single'));
gMap2 = gpuArray(cast(M.cMap2, 'single'));
gpox  = gpuArray(cast(M.poxMap, 'single'));
gscav = gpuArray(cast(M.scavMap,'single'));
gamap = gpuArray(cast(M.axonMap, 'int16'));
gcmap = gpuArray(cast(M.centerMap, 'single'));
%gpixg  = gpuArray(M.pixelMap(:,:,2));
gpixb  = gpuArray(cast(zeros(N),'single'));


fprintf('Running simulation .... \n');

figure();
lowerLimit = 2;

upperLimit = (N)^2-1;
for i=1:simIterations
    %tic
    [gpox, gamap, gpixb, gMap2] = feval( kernel, N, gpox, gamap, gpixb, gMap2, gMap1, gscav,gcmap, 5*M.diffInside, 10*M.diffOutside, lowerLimit, upperLimit ); wait(gd);
    
    %toc
    if mod(i, 200) == 1
        fprintf('Running simulation: [Iter %d of %d]\n', i, simIterations);
        %fprintf("==>> %d\n",i);
          rpixel = gather(gMap2);
%         gpixel = gather(gpix);
         bpixel = gather(gpixb);
         pm = cat(3, rpixel, M.spaceMap, bpixel);
        %                plot_model('pixelmap');
          axis equal, hold on
          imagesc(pm, 'XData', [-1*M.nerve_r M.nerve_r], 'YData', [-1*M.nerve_r, M.nerve_r]);
          axis('on', 'image');
          hold off;
        drawnow
   end
   
   [gpox, gamap, gpixb, gMap1] = feval( kernel, N, gpox, gamap, gpixb, gMap1, gMap2, gscav, gcmap, 5*M.diffInside, 10*M.diffOutside, lowerLimit, upperLimit ); wait(gd);
end
fprintf('Simulation Complete\n');
end