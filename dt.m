function dt
       %x = [-0.5 -0.5 0.5 0.5]';
       %y = [-0.5 0.5 0.5 -0.5]';
       x = (3000*rand(1200,1)-1500);
       y = (3000*rand(1200,1)-1500);
       tri = delaunay(x,y);
 
       % Highlight the first triangle in red
       tri1 = tri(1, [1:end 1]);
       hold on
       plot(x(tri1), y(tri1), 'r')
       
end