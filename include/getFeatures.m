function [features, validPoints] = getFeatures(image)
if length(size(image)) > 2
    baseImage = rgb2gray(image);
else
    baseImage = image;
end
POINTS = detectSURFFeatures(baseImage);
%POINTS = detectHarrisFeatures(baseImage);
%[features, validPoints] = extractFeatures(baseImage, POINTS, 'SURFSize', 128);
[features, validPoints] = extractFeatures(baseImage, POINTS);
