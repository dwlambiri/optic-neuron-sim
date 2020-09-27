z = [
   41.80%
41.40%
41.00%
40.60%
40.70%
40.40%
40.50%
40.40%
40.60%
40.20%
40.20%
39.80%
38.60%
38.80%
38.40%
37.60%
38.40%
38.40%
38.70%
38.20%
38.60%
39.30%
39.30%
37.70%
38.70%
37.80%
39.20%
37.70%
38.70%
37.60%
34.50%
33.90%
34.50%
34.50%
34.10%
36.50%
35.70%
34.60%
27.50%
32.60%
31.60%
35.70%
30.40%
36.00%
24.90%
36.90%
26.50%
34.00%
33.10%
34.90%
21.40%
34.70%
32.30%
35.40%
29.80%
34.00%
20.30%
18.60%
31.70%
22.20% 
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