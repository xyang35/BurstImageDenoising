function [consistentImageSet] = getConsistentImageSet(grayScaleImageSet, homographySet)
consistentImageSet = zeros(size(grayScaleImageSet));

for i = 1 : size(grayScaleImageSet, 3)
    if isempty(homographySet{i})
        consistentImageSet(:,:,i) = grayScaleImageSet(:,:,i);
        continue;
    end
    homographyFlow = homographySet{i};

    consistentImageSet(:,:,i) = backwardTransform(grayScaleImageSet(:,:,i), homographyFlow);
end
