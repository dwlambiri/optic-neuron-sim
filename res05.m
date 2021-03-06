

z = [
    3.60%
2.30%
1.40%
0.90%
0.60%
0.40%
0.30%
0.20%
0.10%
0.10%
6.90%
4.60%
3.30%
2.30%
1.50%
1.10%
0.70%
0.50%
0.30%
0.30%
11.90%
8.10%
5.80%
4.30%
3.10%
2.30%
1.60%
1.20%
0.80%
0.60%
18.80%
13.00%
9.30%
6.90%
5.00%
4.00%
3.00%
2.30%
1.70%
1.30%
29.00%
19.70%
14.10%
10.40%
7.80%
6.00%
4.60%
3.70%
2.90%
2.30%
43.30%
29.00%
20.30%
15.00%
11.40%
8.70%
6.90%
5.40%
4.40%
3.60%
59.40%
41.50%
29.00%
20.80%
15.80%
12.30%
9.60%
7.60%
6.10%
4.90%
72.80%
55.80%
40.30%
29.00%
21.30%
16.40%
13.00%
10.40%
8.30%
6.90%
84.20%
69.00%
54.00%
39.30%
29.00%
21.80%
17.10%
13.70%
11.10%
9.00%
92.10%
80.00%
65.40%
51.10%
38.50%
29.00%
22.20%
17.70%
14.40%
11.90%
96.60%
88.50%
76.00%
62.70%
49.30%
37.90%
29.00%
22.60%
18.30%
15.00%
99.00%
94.20%
85.00%
72.90%
60.50%
47.70%
37.20%
29.00%
23.00%
18.80%
99.70%
97.50%
91.40%
81.60%
69.90%
58.20%
46.50%
36.80%
29.00%
23.30%
99.90%
99.20%
95.70%
88.50%
78.40%
67.00%
55.70%
45.30%
36.30%
29.00%
100.00%
99.70%
98.10%
93.50%
85.50%
75.30%
64.90%
53.90%
44.30%
36.00%
100.00%
99.90%
99.30%
96.60%
91.00%
82.70%
72.80%
62.10%
52.30%
43.30%
100.00%
100.00%
99.80%
98.50%
94.80%
88.40%
80.00%
70.30%
61.00%
51.10%
100.00%
100.00%
99.90%
99.40%
97.30%
92.90%
86.00%
77.20%
68.10%
59.10%
100.00%
100.00%
100.00%
99.80%
98.90%
95.90%
90.60%
83.60%
74.90%
66.00%
100.00%
100.00%
100.00%
99.90%
99.40%
97.90%
94.10%
88.40%
81.10%
72.90%
100.00%
100.00%
100.00%
100.00%
99.80%
99.00%
96.60%
92.40%
86.30%
78.80%
];

xdimsize = 21;
ydimsize = 10;

x = linspace(5, 25,xdimsize);
y = linspace(0.003,0.006, ydimsize);
[X,Y]=meshgrid(x,y);
Z = reshape(z, [ydimsize xdimsize]);
figure; surf(X, Y , Z);
%q = ones(ydimsize, xdimsize)*20;
%hold on; surf(X, Y, q); hold off;
%q1 = ones(ydimsize, xdimsize)*60;
%hold on; surf(X, Y, q1); hold off;
xlabel('death thr zepto/um2'); ylabel('prod zepto/um2/r'); zlabel('perc alive');
title('alive =f(death, prod) @deathExtra = 0');