function fc(file, size, polarity)

onbig = process_image(file, size,[],[]);
[c, r]  = imfindcircles(onbig,[6 60], 'ObjectPolarity', polarity);
hold on; plot(c(:,1)', c(:,2)', 'r.'); hold off;

end