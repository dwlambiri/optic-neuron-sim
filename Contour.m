%% Plots the flooding contour
function upt = Contour(p, t, death_time)
xmin = min(p(1, t)); xmax = max(p(1, t));
ymin = min(p(2, t)); ymax = max(p(2, t));
nt = size(t, 2);
nxy = ceil(sqrt(nt/2)) + 1;
x = linspace(xmin, xmax, nxy);
y = linspace(ymin, ymax, nxy);
[xx, yy] = meshgrid(x, y);
zz=tri2grid(p, t, death_time, x, y);

[~, h] = contourf(xx, yy, max(death_time) - zz);
%colormap(parula)
% caxis([-20,20])

upt = @Upt;
upt(0);
    function Upt(tim)
        h.LevelList = max(death_time) - tim;
    end
end
