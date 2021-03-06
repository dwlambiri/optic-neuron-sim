z = [
    43.10%
42.50%
42.30%
42.10%
41.80%
42.00%
41.00%
40.50%
40.50%
40.60%
40.70%
40.60%
37.80%
37.50%
37.10%
37.00%
37.40%
37.30%
36.50%
35.80%
34.80%
35.80%
35.40%
35.20%
34.50%
34.00%
34.30%
34.50%
33.20%
33.50%
30.70%
30.80%
31.30%
32.80%
30.60%
31.80%
27.00%
30.60%
29.20%
30.30%
28.60%
28.50%
25.70%
27.70%
27.80%
28.70%
27.10%
27.00%
25.10%
25.40%
26.80%
26.70%
26.50%
26.80%
20.40%
23.60%
22.90%
21.70%
23.90%
21.50%

];

xdimsize = 10;
ydimsize = 6;

x = linspace(0.1, 1,xdimsize);
y = linspace(5, 10, ydimsize);
[X,Y]=meshgrid(x,y);
Z = reshape(z, [ydimsize xdimsize]);
figure; surf(X, Y , Z);
xlabel('all diff coef'); ylabel('res pix/um'); zlabel('perc alive');
title('alive =f(diff, res) @thres = 20 prod = 0.005 @deathExtra = 40');