
function propagation_alg_cuda_m01(simIterations, M, itershow, algo, showTox)

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

N = 2*M.opticNerveRadiusPixels;

numPixels = 0;
usedPixels = zeros(1,N*N);

for i=1:N*N
    if M.spaceMap(i) ~= -1
       numPixels = numPixels+1; 
       usedPixels(numPixels) = i -1;
    end
end

ptxFilename =  strcat(filePrefix, '.ptx');
kernel = parallel.gpu.CUDAKernel( ptxFilename, cudaFilename );

kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock,1,1];
kernel.GridSize = [ceil(numPixels/kernel.MaxThreadsPerBlock),1];

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

gMap = gpuArray(cast( cat(3, M.cMap1, M.cMap1), 'single'));
%gMap1 = gMap(:,:,1);
%gMap2 = gMap(:,:,2);
%gMap1 = gpuArray(cast(M.cMap1, 'single'));
gpox  = gpuArray(cast(M.poxMap, 'single'));
gscav = gpuArray(cast(M.scavMap,'single'));
gamap = gpuArray(cast(M.axonMap, 'int16'));
gcmap = gpuArray(cast(M.centerMap, 'single'));
gtmap = gpuArray(cast(M.axonDeathValue, 'single'));
%gpixg  = gpuArray(M.pixelMap(:,:,2));
gpixb  = gpuArray(cast(zeros(N),'uint8'));
gused  = gpuArray(cast(usedPixels(1:numPixels), 'int32'));


fprintf('Running simulation .... \n');

fig = figure();
lowerLimit = 2;

%Q = parallel.pool.DataQueue;
%afterEach(Q,@(data) displayImage(fig, data, M.opticNerveRadiusPixels));

upperLimit = (N)^2-2;
tic
for i=1:simIterations
    %tic
    [gpox, gamap, gpixb, gMap(:,:,2-mod(i,2)), gscav] = feval( kernel, N, numPixels, gused, gpox, gamap, gpixb, gMap(:,:,2- mod(i,2)), gMap(:,:,1+mod(i,2)), gscav,gcmap, M.diffValues(M.diffInside), M.diffValues(M.diffOutside), lowerLimit, upperLimit, M.deathToxThreshold_r, M.extraToxReleaseOnDeath_r, 1-M.scavOutsideAxon_r, algo, gtmap ); 
    
    %toc
   
   %[gpox, gamap, gpixb, gMap(:,:,1), gscav] = feval( kernel, N, gpox, gamap, gpixb, gMap(:,:,1), gMap(:,:,2), gscav, gcmap, M.diffInside, M.diffOutside, lowerLimit, upperLimit, M.deathThr, M.deathRelease, 1 - M.scavOut, algo, gtmap ); 

   if mod(i, itershow) == 1
        fprintf('Running simulation: [Iter %d of %d]\n', i, simIterations);
         %pm = cat(3, gMap1, gMap1, gMap1);
         if showTox
             cla reset; imagesc(gMap(:,:,1)); drawnow('limitrate'); 
         else
             cla reset; imagesc(gpixb); drawnow('limitrate'); 
         end
          
         %displayImage(fig, pm, M.nerve_r);
        %                plot_model('pixelmap');
   end
end
toc
fprintf('Simulation Complete\n');
end