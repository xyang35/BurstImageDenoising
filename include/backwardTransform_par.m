function adjustedImage = backwardTransform_par(ori_image, homographyFlow)
% Paralleled version for backward transformation
% transforming other images to reference image
% loop in  ref image space
% introduce interpolation using neighborhood pixels
parpool(4);

adjustedImage = zeros(size(ori_image));
[rows, cols, ~] = size(homographyFlow);
homographyFlow_r = homographyFlow(:,:,2);
homographyFlow_c = homographyFlow(:,:,1);
parfor r = 1 : rows
    disp(r)
    for c = 1 : cols
        
        r_ = r + homographyFlow_r(r, c);
        c_ = c + homographyFlow_c(r, c);
        
        adjustedImage(r, c, :) = imginterp(ori_image, [r_, c_]);
    end
end

adjustedImage = uint8(adjustedImage);
parpool close

end

function val=imginterp(srcimg, srcpt)

[h, w, d] = size(srcimg);
x1 = min(max(floor(srcpt(1)), 1), h);
y1 = min(max(floor(srcpt(2)), 1), w);
if x1 == h
    x1 = h-1;
    x2 = h;
else
    x2 = x1+1;
end
if y1 == w
    y1 = w-1;
    y2 = w;
else
    y2 = y1+1;
end

[X, Y] = meshgrid(x1:x2, y1:y2);
val = zeros(d, 1);
for depth = 1 : d
    Z = zeros(size(X));
    for i = 1 : numel(Z)
        Z(i) = srcimg(X(i), Y(i), depth);
    end
    val(depth) = interp2(X,Y,Z, srcpt(1), srcpt(2), 'linear', Z(1));
end

end