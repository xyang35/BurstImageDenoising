function homographyFlow = discretizeAndGroupImageHomography(homographySet, nodeRows, nodeCols, imageRows, imageCols)
% Discretize and group image node homographies to form the homography flow
homographyFlow = zeros(imageRows, imageCols, 2);
rowsPerNode = floor(imageRows / nodeRows);
colsPerNode = floor(imageCols / nodeCols);
for r = 1 : nodeRows
    for c = 1 : nodeCols
        homography = homographySet(r, c).homographies;
        rStart = (r - 1) * rowsPerNode + 1;
        rEnd = r * rowsPerNode;
        if r == nodeRows
            rEnd = imageRows;
        end
        cStart = (c - 1) * colsPerNode + 1;
        cEnd = c * colsPerNode;
        if c == nodeCols
            cEnd = imageCols;
        end
        nodeHomographyFLow = discretizeHomography(homography, [rStart rEnd], [cStart cEnd]);
        homographyFlow(rStart : rEnd, cStart : cEnd, :) = nodeHomographyFLow(rStart : rEnd, cStart : cEnd, :);
    end
end
end

function homographyFLow = discretizeHomography(homography, rows, cols)
% The matrix T uses the convention:
% [x y 1] = [u v 1] * T
% where T has the form:
% [a b c;...
%  d e f;...
%  g h i];
%homographyFLow = zeros(rows(2) - rows(1) + 1, cols(2) - cols(1) + 1, 2);
for x = cols(1) : cols(2)
    for y = rows(1) : rows(2)
        pt = [x, y, 1];
        pt_ = pt * homography;
        pt_ = pt_ / pt_(3);
        homographyFLow(y, x, 1) = pt_(1) - pt(1);
        homographyFLow(y, x, 2) = pt_(2) - pt(2);
    end
end
end
