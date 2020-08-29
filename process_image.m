function onibigbw = process_image(imageFile, nerve_r, resize, togrey)
if ~exist( imageFile, 'file')
    fprintf('Cannot load image file %s\n', imageFile);
    return;
else
    oni = imread(imageFile);
    osz = size(oni);
    f1 = nerve_r*2/osz(1);
    if isempty(resize)
        resize = 'bicubic';
    end
    onibig = imresize(oni, f1, resize);
    figure; imshow(onibig);
    
%     [pixelCount, grayLevels] = imhist(onibig);
%     figure; bar(pixelCount);
%     
%     onibig01 = onibig < 40;
%     onibig01 = imfill(onibig01, 'holes');
%     figure; imshow(onibig01);
    
    if isempty(togrey)
        togrey = 'adaptive';
    end
    onibigbw = imbinarize(rgb2gray(onibig),togrey);
    figure; imshow(onibigbw);
end
end