function draw_circles(center, radius, res, doNotFill, allWhite)

x = center(:,1);
y = center(:,2);

%res = min(round(.5*radius), 100);

%if doNotFill 
    for k=1:length(radius)
       if radius(k) < 2
           if allWhite
             plot(x(k),y(k),'w');
           else
             plot(x(k),y(k));  
           end
       else
          ang = (linspace(0,2*pi,max(8,radius(k)/3))); %  
          xp = radius(k)*cos(ang);
          yp = radius(k)*sin(ang);
          xr = x(k)*ones(size(ang))+xp;
          yr = y(k)*ones(size(ang))+yp;
          if allWhite
            if doNotFill
                plot( xr, yr, 'w');
            else
                fill(xr, yr, 'w');
            end
          else 
            if doNotFill
                plot(xr,yr);
            else
                fill(xr, yr, ' ');
            end
          end
       end
    end
% else
%     ang = (linspace(0,2*pi,res)); % 
% 
%     xp = radius*cos(ang);
%     yp = radius*sin(ang);
%     if allWhite
%         fill((x*ones(size(ang))+xp)',(y*ones(size(ang))+yp)','w');
%     else
%        fill((x*ones(size(ang))+xp)',(y*ones(size(ang))+yp)',' '); 
%     end
% end

end