function test_distribution(numPoints)

total = numPoints; 
figure; 
x= skewed_distr(total, 12, 10, 3, 70);

[N, P] = hist(x, 60); 
N = N/total*100; 
hold on; 
yyaxis left; plot(P,N); 
yyaxis right; plot(P,cumsum(N)); 
grid on; 
hold off;