function adjustedImage = backwardTransform(ori_image, homographyFlow)
% Updated version for backward transformation
% transforming other images to reference image
% loop in  ref image space
% introduce interpolation using neighborhood pixels

adjustedImage = zeros(size(ori_image));
[rows, cols, ~] = size(homographyFlow);
for r = 1 : rows
    %disp(r)
    for c = 1 : cols
        
        r_ = round(r + homographyFlow(r, c, 2));
        c_ = round(c + homographyFlow(r, c, 1));
        r_ = min(max(1, r_), rows);
        c_ = min(max(1, c_), cols);
        adjustedImage(r, c, :) = ori_image(r_, c_, :);
    end
end

adjustedImage = uint8(adjustedImage);
end