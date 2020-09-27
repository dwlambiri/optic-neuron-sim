%%
% on death threashold = 2
% insult = 0
% detox extra = 0.001

z = [
85.50%
82.60%
77.60%
71.00%
62.20%
53.50%
46.20%
94.20%
91.60%
84.70%
72.60%
59.90%
49.90%
42.30%
97.10%
95.00%
86.40%
71.40%
60.00%
49.90%
44.40%
98.30%
96.30%
84.90%
68.70%
60.20%
48.70%
39.50%
99.00%
97.20%
84.00%
66.60%
55.90%
44.50%
37.50%
99.40%
97.20%
84.40%
68.60%
54.00%
41.60%
29.90%
99.70%
98.60%
89.50%
65.80%
52.70%
38.70%
13.30%
99.80%
98.50%
89.20%
62.50%
50.80%
35.50%
9.70%
99.80%
98.60%
90.00%
61.80%
48.00%
18.40%
6.00%
99.90%
98.40%
88.00%
58.30%
47.00%
16.70%
3.20%
];

xdimsize = 10;
ydimsize = 7;

x = linspace(0.1, 1,xdimsize)/5;
y = linspace(0,60, ydimsize);
[X,Y]=meshgrid(x,y);
Z = reshape(z, [ydimsize xdimsize]);
figure; surf(X, Y , Z);
xlabel('diffusion coeff'); ylabel('on death SOX zetamol/um2'); zlabel('perc alive');
title('alive =f(diffusion, deathExtra) @prod = 0.0045 zetamol/um2/iter @deathThreash = 20 zeta mol/um2');