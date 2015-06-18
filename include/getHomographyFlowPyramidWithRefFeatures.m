function [homographyFlowPyramid, homographyPyramid] = getHomographyFlowPyramidWithRefFeatures(refPyramid, refFeatures, refPoints, pyramid, FEATURELEVEL)
featImage2 = rgb2gray(pyramid{FEATURELEVEL});

% third, detect feature points and perform matching on lowest layer
% use the points and matching correspondence to compute homography for each
% image node

[features2, validPoints2] = getFeatures(featImage2, FEATURELEVEL);

%POINTS2 = detectSURFFeatures(baseImage2);
% imshow(pyramid1{1});
% hold on
% plot(POINTS1(1 : 100));
%[features2, validPoints2] = extractFeatures(baseImage2, POINTS2);

[indexPairs, ~] = matchFeatures(refFeatures,features2);
disp(size(indexPairs, 1));

matchedPts1 = refPoints(indexPairs(:,1));
matchedPts2 = validPoints2(indexPairs(:,2));
%[~, inlierPts1, inlierPts2] = estimateGeometricTransform(matchedPts1, matchedPts2, 'projective');
%figure; showMatchedFeatures(refPyramid{2},pyramid{2},matchedPts1,matchedPts2);


homographyPyramid = getHomographyPyramid(refPyramid, matchedPts1, matchedPts2, FEATURELEVEL);

% forth, refine homography
%homographyPyramid_old = homographyPyramid;
%homographyPyramid = refineHomographyPyramid(homographyPyramid);
homographyPyramid = refineHomographyPyramid_v2(homographyPyramid);

% fifth, discretize homography to create homography flow
homographyLevel = homographyPyramid{length(homographyPyramid)};
%load('homographyLevel.mat');
imageLevel = pyramid{length(homographyPyramid)};
%homographyFlow = discretizeAndGroupImageHomography_old(homographyLevel, size(homographyLevel, 1), size(homographyLevel, 2), size(imageLevel, 1), size(imageLevel, 2));
homographyFlow = discretizeAndGroupImageHomography(homographyLevel, size(homographyLevel, 1), size(homographyLevel, 2), size(imageLevel, 1), size(imageLevel, 2));

homographyFlowPyramid = cell(1, length(homographyPyramid));
homographyFlowPyramid{length(homographyPyramid)} = homographyFlow;
for level = 1 : length(homographyPyramid) - 1
    imageLevel = pyramid{level};
    homographyFlow = imresize(homographyFlow, [size(imageLevel, 1), size(imageLevel, 2)], 'nearest');
    homographyFlowPyramid{level} = homographyFlow / (2 ^ (length(homographyPyramid) - level));
end
