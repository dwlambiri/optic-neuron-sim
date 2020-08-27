function [n, centers, radii] = make_capillary(x, y, theta, r_big, r_small)
        n =floor(r_big/(2*r_small));
        centers = [linspace(x+r_big*cos(theta)/2, x+r_big * cos(theta), n)' linspace(y+r_big*sin(theta)/2, y+r_big * sin(theta), n)']; 
        radii = linspace(0, r_small, n)';
        
end