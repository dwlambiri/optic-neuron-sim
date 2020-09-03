function test_distribution(numPoints)


epsil = 0.035;

diameter = [0;0.1;0.2;0.3;0.4;0.45;0.5;0.55;0.6;0.65;0.7;0.8;0.9;1;1.1;1.2;1.3;1.4;1.5;1.6;1.7;1.8;1.9;2;2.1;2.2;2.3;2.4;2.5;2.6;2.7;2.8;2.9;3;3.1;3.2;3.3;3.4;3.5;3.6;3.7;3.8;3.9;4;4.1;4.2;4.3;4.4;4.5;4.6;4.7;4.8;4.9;5;5.1;5.2;5.3;5.4;5.5;5.6;5.7;5.8;5.9];
probability = [0;0;0;0.7;1.25;2.95;4;5;6.5;8;9.4;9.6;8.7;7.7;6.7;5.7;4.85;4.1;3.65;3;2.65;2.25;2.1;1.85;1.4;1.3;1.2;1;0.9;0.8;0.65;0.55;0.55;0.45;0.4;0.35;0.3;0.25;0.25;0.2;0.1;0.1;0.1;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil;epsil];

maxidx = find(probability == max(probability));
10*diameter(maxidx)

total = numPoints; 
figure; 
x= skewed_distr(total, 10, 8, 2, 70);

[N, P] = hist(x, 60); 
N = N/total*100; 
hold on; 
yyaxis left; plot(P,N); 
plot(10*diameter, probability);
yyaxis right; plot(P,cumsum(N));

d = 10*[0; diff(diameter)];
dv = d .* probability;
plot(10*diameter, cumsum(dv));
grid on; 
hold off;