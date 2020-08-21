%% Data from the paper

mean_diameter = [1.04 1.19	1.15    1.21    1.28    1.27
                .97  .95    1.01    1.16    1.24    1.26
                .93  .92    1.07    1.16    1.22    1.28
                .94  .93    1.08    1.19    1.18    1.17
                .92  1.02   1.02    1.14    1.16    1.24
                1.01 1.1    1.1     1.22    1.22    1.21];
% figure, imshow(-mean_diameter, [], 'initialmagnification', 6000)

%% Fit using Cubic spline interpolation
[x, y] = meshgrid(-1:(2/5):1);
sf = fit([x(:) y(:)],mean_diameter(:),'cubicinterp');
figure, plot(sf,[x(:) y(:)], mean_diameter(:));

%% Create a fine grid using the fit
grid_siz = 1000;
[x_mesh, y_mesh] = meshgrid(linspace(-1,1, grid_siz));
mean_mesh = single(sf(x_mesh, y_mesh));
mean_mesh_radius = mean_mesh/2;

figure,
plot3(x, y, mean_diameter, 'r.', 'markersize', 20); hold on
surf(x_mesh, y_mesh, mean_mesh, 'linestyle', 'none')

%% Write grid to file

% whos mean_mesh_radius % sse size of data in bytes
fid = fopen('axon_mean_radius.tbl', 'w');
fwrite(fid, grid_siz, 'uint32');
fwrite(fid, mean_mesh_radius, 'single');
fclose(fid);



