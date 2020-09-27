
function propagation_alg_cuda_m01(simIterations, M, itershow, algo)

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
        if M.centerMap(yinsult-M.insultRadius+i, xinsult-M.insultRadius+j) <= 0
            continue;
        end
        M.cMap1(yinsult-M.insultRadius+i, xinsult-M.insultRadius+j) = M.insultAmount;
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

fig = figure();
lowerLimit = 2;

%Q = parallel.pool.DataQueue;
%afterEach(Q,@(data) displayImage(fig, data, M.nerve_r));

upperLimit = (N)^2-2;
tic
for i=1:simIterations/2
    %tic
    [gpox, gamap, gpixb, gMap2, gscav] = feval( kernel, N, gpox, gamap, gpixb, gMap2, gMap1, gscav,gcmap, M.diffInside, M.diffOutside, lowerLimit, upperLimit, M.deathThr, M.deathRelease, 1-M.scavOut, algo ); 
    
    %toc
   
   [gpox, gamap, gpixb, gMap1, gscav] = feval( kernel, N, gpox, gamap, gpixb, gMap1, gMap2, gscav, gcmap, M.diffInside, M.diffOutside, lowerLimit, upperLimit, M.deathThr, M.deathRelease, 1 - M.scavOut, algo  ); 

   if mod(i, itershow) == 1
        fprintf('Running simulation: [Iter %d of %d]\n', 2*i, simIterations);
        %fprintf("==>> %d\n",i);
          %rpixel = gather(gMap2);
          rpixel = zeros(size(M.centerMap));
          bpixel = gather(gpixb);
         %bpixel = gather(gpixb);
         pm = cat(3, rpixel, M.centerMap, bpixel);
         displayImage(fig, pm, M.nerve_r);
        %                plot_model('pixelmap');
   end
end
toc
fprintf('Simulation Complete\n');
end