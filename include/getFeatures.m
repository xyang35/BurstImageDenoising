function [features, validPoints] = getFeatures(image, FEATURELEVEL)
if nargin < 2
    FEATURELEVEL = 1;
end

if length(size(image)) > 2
    grayImage = rgb2gray(image);
else
    grayImage = image;
end

POINTS = detectSURFFeatures(grayImage);

% filter feature points smaller than scale
%scale = 1.6 * FEATURELEVEL;
%POINTS = POINTS(find(POINTS.Scale > scale));

%POINTS = detectSURFFeatures(grayImage);
%POINTS = detectHarrisFeatures(baseImage);
%[features, validPoints] = extractFeatures(baseImage, POINTS, 'SURFSize', 128);
[features, validPoints] = extractFeatures(grayImage, POINTS);
