function test_fit(howmany, w1, w2, w3)

[b, p, c] =pdf_pan_paper(1);


x = zeros(1, howmany);

for i=1:howmany
    diam = howmany/600;
    x(i) = p(diam, w1, w2, w3);
end

size(x) 
[N, P] = hist(x, 60); 
N = N/howmany*100; 
figure;
hold on; 
yyaxis left; plot(P,N); 
yyaxis right; plot(P,cumsum(N)); 
grid on; 
hold off;

end