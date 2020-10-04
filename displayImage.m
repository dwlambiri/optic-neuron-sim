function displayImage(fig, pm, nerve_r)

%figure(fig);
axis equal, hold on
cla reset;
imshow(pm, 'XData', [-1*nerve_r nerve_r], 'YData', [-1*nerve_r, nerve_r]);
axis('on', 'image');
hold off;
drawnow('limitrate');         
end