
function propagation_alg_cuda_m3d(noPlanes, simIterations, M, itershow, algo, showTox, injurySize)

gd = gpuDevice();
reset(gd);
wait(gd);

% Propagation based on closest points and length of connections

%neuron_speed_formula = @(R) 2./R; % 2/R, For  2/R^2  write this code   ->>  2./R.^2
% Load the kernel
filePrefix = 'propagation_algo_cuda3d';
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

noExtra  = 2;
d3 = cat(3, M.cMap1, M.cMap1);
for i=1:noPlanes
    d3 = cat(3,d3,M.cMap1);
end
gMap = gpuArray(cast( d3, 'single'));
clear d3;
%gMap1 = gMap(:,:,1);
%gMap2 = gMap(:,:,2);
%gMap1 = gpuArray(cast(M.cMap1, 'single'));
%gMap2 = gpuArray(cast(M.cMap2, 'single'));
gpox  = gpuArray(cast(M.poxMap, 'single'));
gscav = gpuArray(cast(M.scavMap,'single'));
gamap = gpuArray(cast(M.axonMap, 'int16'));
gcmap = gpuArray(cast(M.centerMap, 'single'));
gtmap = gpuArray(cast(M.axonDeathValue, 'single'));
%gpixg  = gpuArray(M.pixelMap(:,:,2));
gpixb  = gpuArray(cast(zeros(N),'uint8'));


fprintf('Running simulation .... \n');

figure();
lowerLimit = 2;

%Q = parallel.pool.DataQueue;
%afterEach(Q,@(data) displayImage(fig, data, M.nerve_r));

upperLimit = (N)^2-2;
totalPlanes = noPlanes + noExtra;
tic

head = 2;
% for i=0:simIterations
%     %tic
%     fprintf("LEVEL %d [%d, %d, %d, %d]\n", 1, 1+mod(head-2, totalPlanes), 1+mod(head, totalPlanes), 1+mod(head, totalPlanes), 1+mod(head+1, totalPlanes));    
%     for j = 1: noPlanes-2
%         fprintf("LEVEL %d [%d, %d, %d, %d %d]\n", j, 1+mod(head-2+j, totalPlanes), 1+mod(head-2+j, totalPlanes), 1+mod(head-1+j, totalPlanes), 1+mod(head+j, totalPlanes), 1+mod(head+j+1, totalPlanes));
%     end
%     fprintf("LEVEL %d [%d, %d, %d, %d %d]\n\n\n", noPlanes, 1+mod(head-2+noPlanes-1, totalPlanes), 1+mod(head-2+noPlanes-1, totalPlanes), 1+mod(head-1+noPlanes-1, totalPlanes), 1+mod(head+noPlanes-1, totalPlanes), 1+mod(head+noPlanes-1, totalPlanes));
% 
%     head = mod(head-2, totalPlanes);
% end


for i=0:simIterations
    %tic
        %fprintf("LEVEL %d [%d, %d, %d, %d]\n", 1, 1+mod(head-2, totalPlanes), 1+mod(head, totalPlanes), 1+mod(head, totalPlanes), 1+mod(head+1, totalPlanes));    

    [gpox, gamap, gpixb, gMap(:,:,1+mod(head-2, totalPlanes)), gscav] = feval( kernel, N, gpox, gamap, gpixb, gMap(:,:,1+mod(head-2, totalPlanes)), gMap(:,:,1+mod(head, totalPlanes)), gMap(:,:,1+mod(head, totalPlanes)), gMap(:,:,1+mod(head+1, totalPlanes)), gscav,gcmap, M.diffInside, M.diffOutside, lowerLimit, upperLimit, M.deathThr, M.deathRelease, 1-M.scavOut, algo, gtmap, 1, 0, 0 ); 
    for j = 1: noPlanes-2
        insult = 0;
        if j > 4
            insult = 1;
        end
        %fprintf("LEVEL %d [%d, %d, %d, %d %d]\n", j, 1+mod(head-2+j, totalPlanes), 1+mod(head-2+j, totalPlanes), 1+mod(head-1+j, totalPlanes), 1+mod(head+j, totalPlanes), 1+mod(head+j+1, totalPlanes));
        [gpox, gamap, gpixb, gMap(:,:, 1+mod(head-2+j, totalPlanes)), gscav] = feval( kernel, N, gpox, gamap, gpixb, gMap(:,:,1+mod(head-2+j, totalPlanes)), gMap(:,:, 1+mod(head-1+j, totalPlanes)), gMap(:,:, 1+mod(head+j, totalPlanes)), gMap(:,:, 1+mod(head+j+1, totalPlanes)), gscav,gcmap, M.diffInside, M.diffOutside, lowerLimit, upperLimit, M.deathThr, M.deathRelease, 1-M.scavOut, algo, gtmap, 0, 0, insult); 
    end
    [gpox, gamap, gpixb, gMap(:,:, 1+mod(head-2+noPlanes-1, totalPlanes)), gscav] = feval( kernel, N, gpox, gamap, gpixb, gMap(:,:,1+mod(head-2+noPlanes-1, totalPlanes)), gMap(:,:, 1+mod(head-1+noPlanes-1, totalPlanes)), gMap(:,:, 1+mod(head+noPlanes-1, totalPlanes)), gMap(:,:, 1+mod(head+noPlanes-1, totalPlanes)), gscav,gcmap, M.diffInside, M.diffOutside, lowerLimit, upperLimit, M.deathThr, M.deathRelease, 1-M.scavOut, algo, gtmap, 0, 1, 0 );

    %fprintf("LEVEL %d [%d, %d, %d, %d %d]\n\n\n", noPlanes, 1+mod(head-2+noPlanes-1, totalPlanes), 1+mod(head-2+noPlanes-1, totalPlanes), 1+mod(head-1+noPlanes-1, totalPlanes), 1+mod(head+noPlanes-1, totalPlanes), 1+mod(head+noPlanes-1, totalPlanes));

    head = mod(head-2, totalPlanes);
    
    %toc
   
   if mod(i, itershow) == 1
        fprintf('Running simulation: [Iter %d of %d]\n', i, simIterations);
         %pm = cat(3, gMap1, gMap1, gMap1);
         if showTox
             cla reset; imagesc(gMap(:,:,noPlanes - mod(noPlanes+i, noPlanes))); drawnow('limitrate'); 
         else
             cla reset; imagesc(gpixb); drawnow('limitrate'); 
         end
          
         %displayImage(fig, pm, M.nerve_r);
        %                plot_model('pixelmap');
   end
end
toc
fprintf('Simulation Complete. HEAD = %d\n', head);
toxMap = gather(gMap);
imshow3D(toxMap);
drawnow('limitrate'); 

% for i = 1: noPlanes+noExtra
%    sum(toxMap(:,:,i), 'all') 
% end
end