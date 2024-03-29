%%
% all diff rates = 0.1
% on death extra = 6
% on death threashold = 2
% insult = 0
% detox extra = 0.001
z= [
 100.00%
100.00%
100.00%
100.00%
100.00%
100.00%
100.00%
99.10%
100.00%
100.00%
100.00%
100.00%
100.00%
100.00%
80.50%
99.40%
100.00%
100.00%
100.00%
100.00%
100.00%
53.40%
90.70%
99.60%
100.00%
100.00%
100.00%
100.00%
35.40%
64.80%
95.10%
99.70%
100.00%
100.00%
100.00%
22.70%
45.30%
74.50%
96.00%
99.80%
100.00%
100.00%
16.00%
30.70%
54.20%
82.50%
97.00%
99.80%
100.00%
11.70%
22.00%
40.20%
63.00%
87.20%
97.50%
99.70%
8.90%
16.60%
29.30%
47.50%
70.30%
89.30%
97.70%
7.00%
12.80%
21.40%
35.60%
55.20%
77.60%
91.70%
5.70%
10.20%
16.90%
27.80%
42.50%
63.30%
81.90%
4.90%
8.50%
13.70%
21.20%
33.30%
49.30%
68.90%
4.30%
7.10%
11.40%
17.20%
26.70%
38.90%
56.00%
3.80%
6.00%
9.60%
14.30%
21.10%
31.50%
44.60%
3.40%
5.30%
8.30%
12.50%
17.70%
25.90%
36.40%
3.00%
4.70%
7.20%
10.80%
14.90%
21.10%
30.20%
2.50%
4.10%
6.20%
9.00%
12.80%
17.30%
24.20%
2.20%
3.70%
5.30%
7.70%
11.00%
14.60%
19.40%
1.90%
3.20%
4.60%
6.80%
9.60%
12.80%
16.40%
1.50%
2.70%
4.10%
5.90%
8.20%
11.30%
14.50%
];


xdimsize = 20;
ydimsize = 7;


x = linspace(0.00233,0.009, xdimsize);
y = linspace(0.00075, 0.0015,ydimsize);

[X,Y]=meshgrid(x,y);
figure; surf(X, Y ,reshape(z, [ydimsize xdimsize]));
xlabel('tox prod'); ylabel('scav'); zlabel('perc alive');
title('alive =f(tox, scav) @deathExtra = 6  @deathThreash = 2');