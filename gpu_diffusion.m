function [c, p, r,g,b] = gpu_diffusion(cMap1, cMap1ym , cMap1yp, cMap1xm, cMap1xp, pox, scav, smap, di)

if smap < 0
    r = 0;
    g = 0;
    b = 0;
    c = 0;
    p = 0;
    return;
end

r= 0;
if smap > 0
    if pox > 0
        g = 1;
        b = 0;
    else
        g = 0;
        b = 1;
    end
    p = pox;
else
    p = 0;
    g = 0;
    b = 0;
end

c = (cMap1+ pox + di*(cMap1ym- cMap1) +  di*(cMap1yp- cMap1) +  di*(cMap1xm-cMap1) + di*(cMap1xp- cMap1))*scav ;
if c > 1
    r = 1;
end
if smap > 0 && pox > 0 && c > 22
    g = 0;
    b = 1;
    p = 0;
    c = 10000;
    
end
end