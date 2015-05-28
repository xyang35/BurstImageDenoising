function pyramid = getPyramids(image)
% Get Gaussian image pyramids from the input image
minRows = 400;
minCols = 400;
[rows, cols, ~] = size(image);
layerNum = 1;
%while rows > minRows && cols > minCols
while rows > minRows || cols > minCols    % fix bug 1
    layerNum = layerNum + 1;
    rows = ceil(rows / 2);
    cols = ceil(cols / 2);
end
pyramid = cell(1, layerNum);
current = image;
pyramid{layerNum} = current;
for level = layerNum - 1 : -1 : 1
    current = impyramid(current, 'reduce');
    pyramid{level} = current;
end
