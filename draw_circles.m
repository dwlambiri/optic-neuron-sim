function draw_circles(center, radius, res, empty, floodThreshold, varargin)

x = center(:,1);
y = center(:,2);

ang = (linspace(0,2*pi,res)); % min(round(.5*radius), 15)
xp = radius*cos(ang);
yp = radius*sin(ang);
if empty 
    h= plot((x*ones(size(ang))+xp)',(y*ones(size(ang))+yp)', varargin{:});
else
    h= fill((x*ones(size(ang))+xp)',(y*ones(size(ang))+yp)',' ', varargin{:});
end

for k = 1:length(h)
    if radius(k) > floodThreshold
        h(k).ZData = ones(size(ang)); 
    end
end

end