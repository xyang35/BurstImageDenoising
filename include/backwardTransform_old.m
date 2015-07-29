function adjustedImage = backwardTransform_old(ori_image, homographyFlow)
% this is a backward transformation, transforming other images to
% reference image
adjustedImage = ori_image;
[rows, cols, ~] = size(homographyFlow);
for r = 1 : rows
    for c = 1 : cols
        r_ = floor(r - homographyFlow(r, c, 2));
        c_ = floor(c - homographyFlow(r, c, 1));
        r_ = min(max(1, r_), rows);
        c_ = min(max(1, c_), cols);
        adjustedImage(r_, c_, :) = ori_image(r, c, :);
    end
end
adjustedImage = uint8(adjustedImage);
end