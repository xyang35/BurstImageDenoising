function homographyPyramid = getHomographyPyramid(pyramid, matchedPoints1, matchedPoints2)
% Compute homography pyramid from matched feature points
featuredPyramid1 = getFeaturedPyramid(pyramid, matchedPoints1);
featuredPyramid2 = getFeaturedPyramidWithRef(pyramid, matchedPoints2, featuredPyramid1);
[~, layerNum] = size(pyramid);
homographyPyramid = cell(1, layerNum);
for level = 1:layerNum
    nodes1 = featuredPyramid1{level};
    nodes2 = featuredPyramid2{level};
    [rows, cols] = size(nodes1);
    homographies = cell(rows, cols);
    for r = 1:rows
        for c = 1:cols
            if length(nodes1{r,c}.pts) < 8 && level > 1% minimum number of points
                r_upper = ceil(r / 2);
                c_upper = ceil(c / 2);
                temp = homographyPyramid{level - 1};
                homographies(r,c) = {temp(r_upper, c_upper).homographies};
            else
                [tform,inlierPts1,inlierPts2] = estimateGeometricTransform(nodes1{r,c}.pts, nodes2{r,c}.pts, 'projective');
                homographies(r,c) = {tform.T};
            end
        end
    end
    pointNum = length(inlierPts1);
    current = struct('homographies', homographies, 'pointNumber', pointNum);
    homographyPyramid(level) = {current};
end
end

function featuredPyramid = getFeaturedPyramidWithRef(pyramid, pts, refFeaturedPyramid)
[~, layerNum] = size(pyramid);
featuredPyramid = cell(1, layerNum);
for level = 1:layerNum
    rows = 2 ^ (level - 1);
    cols = rows;
    multi = rows;
    nodes = cell(rows, cols);
    refNodes = refFeaturedPyramid{level};
    % node: range, featurePoints, homography
    nodePixelNumVert = floor(size(pyramid{level}, 1) / multi);
    nodePixelNumHori = floor(size(pyramid{level}, 2) / multi);
    newPts = pts;
    newPts.Scale = newPts.Scale * multi;
    newPts.Location = newPts.Location * multi;
    for r = 1 : rows
        for c = 1 : cols
            range = [1 + (c - 1) * nodePixelNumHori, c * nodePixelNumHori; 1 + (r - 1) * nodePixelNumVert, r * nodePixelNumVert];
            refNode = refNodes{r,c};
            node = ImageNode(range, newPts(refNode.inds), refNode.inds);
            nodes(r,c) = {node};
        end
    end
    featuredPyramid(level) = {nodes};
end
end

function featuredPyramid = getFeaturedPyramid(pyramid, pts)
[~, layerNum] = size(pyramid);
featuredPyramid = cell(1, layerNum);
for level = 1:layerNum
    rows = 2 ^ (level - 1);
    cols = rows;
    multi = rows;
    nodes = cell(rows, cols);
    % node: range, featurePoints, homography
    nodePixelNumVert = floor(size(pyramid{level}, 1) / multi);
    nodePixelNumHori = floor(size(pyramid{level}, 2) / multi);
    newPts = pts;
    newPts.Scale = newPts.Scale * multi;
    newPts.Location = newPts.Location * multi;
    for r = 1 : rows
        for c = 1 : cols
            range = [1 + (c - 1) * nodePixelNumHori, c * nodePixelNumHori; 1 + (r - 1) * nodePixelNumVert, r * nodePixelNumVert];
            inds = [];
            for i = 1 : length(newPts)
                if inRange(range, newPts(i).Location)
                    inds = [inds, i];
                end
            end
            nodePts = newPts(inds);
            node = ImageNode(range, nodePts, inds);
            nodes(r, c) = {node};
        end
    end
    featuredPyramid(level) = {nodes};
end
end

function flag = inRange(range, pt)
flag = false;
if pt(1) >= range(1,1) && pt(1) < range(1,2) && pt(2) >= range(2,1) && pt(2) < range(2,2)
    flag = true;
end
end