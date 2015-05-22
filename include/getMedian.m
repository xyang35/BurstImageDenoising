function medianImage = getMedian(grayScaleImageSet, homographySet)
homographySet = floor(homographySet);
consistentImageSet = zeros(size(grayScaleImageSet));
rows = size(grayScaleImageSet, 1);
cols = size(grayScaleImageSet, 2);
for i = 1 : size(grayScaleImageSet, 3)
    for r = 1 : rows
        for c = 1 : cols
            ri = min(max(1, r + homographySet(r,c,2,i)), rows);
            ci = min(max(1, c + homographySet(r,c,1,i)), cols);
            consistentImageSet(r,c,i) = grayScaleImageSet(ri,ci,i);
        end
    end
end
medianImage = median(consistentImageSet, 3);