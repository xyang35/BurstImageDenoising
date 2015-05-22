function homographyFlowPyramid = getHomographyFlowPyramid(refPyramid, pyramid)
baseImage1 = rgb2gray(refPyramid{1});
baseImage2 = rgb2gray(pyramid{1});

% third, detect feature points and perform matching on lowest layer
% use the points and matching correspondence to compute homography for each
% image node

[features1, validPoints1] = getFeatures(baseImage1);
[features2, validPoints2] = getFeatures(baseImage2);

%POINTS1 = detectSURFFeatures(baseImage1);
%POINTS2 = detectSURFFeatures(baseImage2);
% imshow(pyramid1{1});
% hold on
% plot(POINTS1(1 : 100));
%[features1, validPoints1] = extractFeatures(baseImage1, POINTS1);
%[features2, validPoints2] = extractFeatures(baseImage2, POINTS2);

[indexPairs, ~] = matchFeatures(features1,features2);
matchedPts1 = validPoints1(indexPairs(:,1));
matchedPts2 = validPoints2(indexPairs(:,2));
%[~, inlierPts1, inlierPts2] = estimateGeometricTransform(matchedPts1, matchedPts2, 'projective');

homographyPyramid = getHomographyPyramid(refPyramid, matchedPts1, matchedPts2);

% forth, refine homography
homographyPyramid = refineHomographyPyramid(homographyPyramid);

% fifth, discretize homography to create homography flow
homographyLevel = homographyPyramid{length(homographyPyramid)};
imageLevel = pyramid{length(homographyPyramid)};
homographyFlow = discretizeAndGroupImageHomography(...
    homographyLevel, size(homographyLevel, 1), size(homographyLevel, 2), size(imageLevel, 1), size(imageLevel, 2));

homographyFlowPyramid = cell(1, length(homographyPyramid));
homographyFlowPyramid{length(homographyPyramid)} = homographyFlow;
for level = 1 : length(homographyPyramid) - 1
    imageLevel = pyramid{level};
    homographyFlow = imresize(homographyFlow, [size(imageLevel, 1), size(imageLevel, 2)], 'nearest');
    homographyFlowPyramid{level} = homographyFlow;
end
