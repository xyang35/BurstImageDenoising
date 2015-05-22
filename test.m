%close all; clear all;
addpath(genpath('./include'));

%% Initialization

% import burst of images
imageNum = 10;
imageSet = cell(1, imageNum);

base_dir = '../burstimages_v1/';
name = 'Bookshelf_2';
path = [base_dir, name];
result_path = ['/localdisk/xyang/PS_data/', name];

ratio = 1;
for i = 1 : imageNum
    image_dir = fullfile(path, [num2str(i - 1), '.jpg']);
    imageSet{i} = imresize(imread(image_dir), ratio);
end
% set the reference image to be the 5th one
ref = 5;
refImage = imageSet{ref}; % get reference image
refPyramid = getPyramids(refImage); % get pyramid for reference image
[refFeatures, refPoints] = getFeatures(refPyramid{1}); % get features and valid points of reference image

%% Homography flows

% homography flows for all base images (base: coarsest level of the pyramid)
baseHomographySet = zeros(size(refPyramid{1}, 1), size(refPyramid{1}, 2), 2, length(imageSet));
% all base images
baseImageSet = cell(1, length(imageSet));
homographyFlowPyramidSet = cell(1, length(refPyramid));
for j = 1 : length(refPyramid)
    homographyFlowPyramidSet{j} = cell(1,length(imageSet));
end
for i = 1 : length(imageSet)
    if i == ref
        baseImageSet{i} = refPyramid{1};
        continue;
    end
    pyramid = getPyramids(imageSet{i});
    homographyFlowPyramid = getHomographyFlowPyramidWithRefFeatures(refPyramid, refFeatures, refPoints, pyramid);
    for j = 1 : length(homographyFlowPyramid)
        homographyFlowPyramidSet{j}{i} = homographyFlowPyramid{j};
    end
    baseHomographySet(:,:,:,i) = homographyFlowPyramid{1};
    baseImageSet{i} = pyramid{1};
    disp(['homography ', num2str(i), ' complete']);
    %return
end
%save([result_path, '_h_norefine.mat'], 'homographyFlowPyramidSet');

%% Consistent pixel

% take grayscale base images to compute pixel difference
baseGrayScaleImageSet = zeros(size(refPyramid{1}, 1), size(refPyramid{1}, 2), length(imageSet));
for i = 1 : length(imageSet)
    baseGrayScaleImageSet(:, :, i) = rgb2gray(baseImageSet{i});
end
% compute median image of all base images (for median value based consistent pixel selection)
[baseConsistentImageSet, ConsistentPixelMap] = getConsistentImageSet(baseGrayScaleImageSet, homographyFlowPyramidSet{1});
baseMedianImage = median(baseConsistentImageSet, 3);
baseIntegralMedImage = integralImage(baseMedianImage); % integral image of the median image
% integral images of all base images
baseIntegralImageSet = zeros(size(refPyramid{1}, 1) + 1, size(refPyramid{1}, 2) + 1, length(imageSet));
for i = 1 : length(imageSet)
    baseIntegralImageSet(:, :, i) = integralImage(baseGrayScaleImageSet(:, :, i));
end
% integral images of all consistent images
ConsistentIntegralImageSet = zeros(size(baseConsistentImageSet, 1) + 1, size(baseConsistentImageSet, 2) + 1, size(baseConsistentImageSet, 3));
for i = 1 : size(ConsistentIntegralImageSet,3)
    ConsistentIntegralImageSet(:, :, i) = integralImage(baseConsistentImageSet(:, :, i));
end

tau = 10; % threshold for selecting consistent pixels
% set of consistent pixel indexes: reference based and median based
baseRefConsistentPixelMap = zeros(size(baseImageSet{1}, 1), size(baseImageSet{1}, 2), length(baseImageSet));
baseMedConsistentPixelMap = zeros(size(baseImageSet{1}, 1), size(baseImageSet{1}, 2), length(baseImageSet));
% set of potential consistent pixel indexes
halfWidth = 2; halfHeight = 2;    % for patch implementation
rows = size(baseRefConsistentPixelMap, 1);
cols = size(baseRefConsistentPixelMap, 2);

% rs = 1 : 10 : rows; cs = 1 : 10 : cols;
% [x, y] = meshgrid(cs, rs);
% quiver(x,y,baseHomographySet(rs,cs,1,i),baseHomographySet(rs,cs,2,i));


for r = 1 : rows
    for c = 1 : cols
        % first get the pixel from the reference image
        sR = max(1, r - halfHeight);
        sC = max(1, c - halfWidth);
        eR = min(rows, r + halfHeight);
        eC = min(cols, c + halfWidth);
        pixNum = (eR - sR + 1) * (eC - sC + 1);
        refPix = baseIntegralImageSet(eR+1,eC+1,ref) - baseIntegralImageSet(eR+1,sC,ref) - baseIntegralImageSet(sR,eC+1,ref) + baseIntegralImageSet(sR,sC,ref);
        refPix = refPix / pixNum;
        medPix = baseIntegralMedImage(eR+1,eC+1) - baseIntegralMedImage(eR+1,sC) - baseIntegralMedImage(sR,eC+1) + baseIntegralMedImage(sR,sC);
        medPix = medPix / pixNum;

        % reference-based
        baseRefConsistentPixelMap(r,c,ref) = 1;
        % from reference to left
        for i = ref - 1 : -1 : 1
            % make use of ConsistenImage to get consistent pixels
            iPix = ConsistentIntegralImageSet(eR+1,eC+1,i) - ConsistentIntegralImageSet(eR+1,sC,i) - ConsistentIntegralImageSet(sR,eC+1,i) + ConsistentIntegralImageSet(sR,sC,i);
            iPix = iPix / pixNum;
            if abs(refPix - iPix) < tau
                baseRefConsistentPixelMap(r,c,i) = 1;
            else
                break;
            end
        end
        % from reference to right
        for i = ref + 1 : length(imageSet)
            iPix = ConsistentIntegralImageSet(eR+1,eC+1,i) - ConsistentIntegralImageSet(eR+1,sC,i) - ConsistentIntegralImageSet(sR,eC+1,i) + ConsistentIntegralImageSet(sR,sC,i);
            iPix = iPix / pixNum;
            if abs(refPix - iPix) < tau
                baseRefConsistentPixelMap(r,c,i) = 1;
            else
                break;
            end
        end

        % median-based
        for i = 1 : length(imageSet)
            iPix = ConsistentIntegralImageSet(eR+1,eC+1,i) - ConsistentIntegralImageSet(eR+1,sC,i) - ConsistentIntegralImageSet(sR,eC+1,i) + ConsistentIntegralImageSet(sR,sC,i);
            iPix = iPix / pixNum;
            if abs(medPix - iPix) < tau
                baseMedConsistentPixelMap(r,c,i) = 1;
            end
        end
    end
end

% combine median consistent pixels and reference consistent pixels
baseConsistentPixelMap = zeros(size(baseRefConsistentPixelMap));
reliableNumber = floor(imageNum / 2);
consistentPixelNumMap = sum(baseMedConsistentPixelMap, 3);
consistentPixelNumMap = consistentPixelNumMap > reliableNumber;
% perform majority filter
consistentPixelNumMap = bwmorph(consistentPixelNumMap, 'majority');

for r = 1 : rows
    for c = 1 : cols
        % case 1: union of the two
        if baseMedConsistentPixelMap(r,c,ref) == 1
            baseConsistentPixelMap(r,c,:) = baseRefConsistentPixelMap(r,c,:) | baseMedConsistentPixelMap(r,c,:);
        % case 2: judge if median based is reliable
        elseif consistentPixelNumMap(r,c) == 1
            % median based result is reliable
            baseConsistentPixelMap(r,c,:) = baseMedConsistentPixelMap(r,c,:);
        else
            % median based result is not reliable
            baseConsistentPixelMap(r,c,:) = baseRefConsistentPixelMap(r,c,:);
        end
    end
end

sumBaseConsistentPixelMap = sum(baseConsistentPixelMap, 3);
[rs, cs] = find(sumBaseConsistentPixelMap == 0);
baseConsistentPixelMap(rs, cs, ref) = 1;

% for i = 1 : size(baseConsistentPixelMap,3)
%     baseConsistentPixelMap(:, :, i) = bwmorph(baseConsistentPixelMap(:, :, i), 'majority');
% end

%% fusion stage

% first estimate noise
fineRefImage = refPyramid{length(refPyramid)};
fineGrayScaleRefImage = double(rgb2gray(fineRefImage));
nontextureMap = imresize(edge(baseMedianImage), size(fineGrayScaleRefImage), 'near');
inds = find(nontextureMap == 0);
fineMedianImage = imresize(baseMedianImage, size(fineGrayScaleRefImage), 'bilinear');
sigma2 = computeSigma2FromDiffVector(fineGrayScaleRefImage(inds) - fineMedianImage(inds));

% second perform temporal fusion
baseConsistentPixelMapR = ConsistentPixelMap(:,:,:,1) .* baseConsistentPixelMap;
baseConsistentPixelMapC = ConsistentPixelMap(:,:,:,2) .* baseConsistentPixelMap;
sigmat2MapSet = cell(1,length(refPyramid));
for level = 1 : length(refPyramid)
    rows = size(refPyramid{level}, 1);
    cols = size(refPyramid{level}, 2);
    % reuse the consistent pixel map for all levels
    if level > 1
        levelConsistentPixelMapR = imresize(baseConsistentPixelMapR, [rows, cols], 'near');
        levelConsistentPixelMapC = imresize(baseConsistentPixelMapC, [rows, cols], 'near');
    else
        levelConsistentPixelMapR = baseConsistentPixelMapR;
        levelConsistentPixelMapC = baseConsistentPixelMapC;
    end
    % get the set of all frames at this level
    levelImageSet = [];
    for i = 1 : length(imageSet)
        if i == ref
            levelImageSet = cat(4, levelImageSet, refPyramid{level});
            continue;
        else
            if level < length(refPyramid)
                ithImage = imresize(imageSet{i},[rows,cols], 'bilinear');
            else
                ithImage = imageSet{i};
            end
            ithImage = backwardTransform(ithImage, homographyFlowPyramidSet{level}{i});%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%
            levelImageSet = cat(4, levelImageSet, ithImage);
        end
    end
    % use consistent image to compute mean value and variance
    levelGrayScaleImageSet = zeros(size(levelImageSet,1), size(levelImageSet,2), size(levelImageSet,4));
    for i = 1 : length(imageSet)
        levelGrayScaleImageSet(:,:,i) = rgb2gray(levelImageSet(:,:,:,i));
    end
    % get consistent image set
    levelConsistentImageSet = levelGrayScaleImageSet; %getConsistentImageSet(levelGrayScaleImageSet, homographyFlowPyramidSet{level});
    levelConsistentPixelMap = levelConsistentPixelMapR > 0;
    levelConsistentImageSet = levelConsistentImageSet .* levelConsistentPixelMap;
    meanImage = sum(levelConsistentImageSet, 3) ./ sum(double(levelConsistentPixelMap), 3);
    % sigma_t^2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigmat2MapSet{level} = sum((levelConsistentImageSet - repmat(meanImage, [1 1 size(levelConsistentImageSet, 3)])) .^ 2, 3) ./ sum(levelConsistentPixelMap, 3);
    % sigma_c^2
    sigmac2Map = sigmat2MapSet{level} - sigma2;
    sigmac2Map = (sigmac2Map > 0) .* sigmac2Map;
    sigmac2Map = sigmac2Map ./ (sigmac2Map + sigma2);
    sigmac2Map = repmat(sigmac2Map, [1 1 3]);
    levelConsistentPixelMap = reshape(levelConsistentPixelMap,...
        [size(levelConsistentPixelMap,1), size(levelConsistentPixelMap,2), 1, size(levelConsistentPixelMap,3)]);
    levelConsistentPixelMap = repmat(levelConsistentPixelMap, [1,1,3,1]);
    meanImage = sum(double(levelImageSet) .* double(levelConsistentPixelMap), 4) ./ sum(double(levelConsistentPixelMap), 4);
    refPyramid{level} = meanImage + sigmac2Map .* (double(refPyramid{level}) - meanImage);
end

% compute sigma for the finest level reference image
differenceImage = [];
rs = [1 2 2 3];
cs = [2 1 3 2];
for i = 1 : 4
    h = zeros(3,3);
    h(2,2) = 1;
    h(rs(i), cs(i)) = -1;
    differenceImage = cat(3, differenceImage, abs(conv2(fineGrayScaleRefImage, h, 'same')));
end
differenceImage = max(differenceImage, [], 3);
differenceImage = 1 ./ (1 + exp(-5 * differenceImage / (sqrt(sigma2) - 3)));
differenceImage = repmat(differenceImage, [1 1 3]);
for level = 2 : length(refPyramid)
    levelSpatiallyFilteredImage = refPyramid{level};
    [rows, cols, ~] = size(levelSpatiallyFilteredImage);
    halfWidth = 2; halfHeight = 2;
    levelGrayScaleImage = double(rgb2gray(uint8(levelSpatiallyFilteredImage)));
    levelGrayScaleGradientImage = zeros(rows, cols, 2);
    hy = [1 2 1; 0 0 0; -1 -2 -1];
    hx = [-1 0 1; -2 0 2; -1 0 1];
    levelGrayScaleGradientImage(:,:,1) = conv2(levelGrayScaleImage, hx, 'same'); % vertical
    levelGrayScaleGradientImage(:,:,2) = conv2(levelGrayScaleImage, hy, 'same'); % horizontal
    levelGrayScaleGradientImage = atan2(levelGrayScaleGradientImage(:,:,2), levelGrayScaleGradientImage(:,:,1));
    levelGrayScaleGradientImage = levelGrayScaleGradientImage / pi * 180;
    for r = 1 + halfHeight : rows - halfHeight
        for c = 1 + halfWidth : cols - halfWidth
            if (levelGrayScaleGradientImage(r,c) >= -22.5 && levelGrayScaleGradientImage(r,c) < 22.5)...
                    || levelGrayScaleGradientImage(r,c) > 157.5 || levelGrayScaleGradientImage(r,c) + 180 < 22.5
                % vertical
                spatialConsistentPixels = cat(4, levelSpatiallyFilteredImage(r-2,c,:), levelSpatiallyFilteredImage(r-1,c,:),...
                    levelSpatiallyFilteredImage(r,c,:), levelSpatiallyFilteredImage(r+1,c,:), levelSpatiallyFilteredImage(r+2,c,:));
%                 levelGrayScaleGradientImage(r,c) = 1;%%
            end
            if (levelGrayScaleGradientImage(r,c) >= 22.5 && levelGrayScaleGradientImage(r,c) <= 67.5)...
                    || (levelGrayScaleGradientImage(r,c) + 180 >= 22.5 && levelGrayScaleGradientImage(r,c) + 180 < 67.5)
                % main diagonal
                spatialConsistentPixels = cat(4, levelSpatiallyFilteredImage(r-2,c-2,:), levelSpatiallyFilteredImage(r-1,c-1,:),...
                    levelSpatiallyFilteredImage(r,c,:), levelSpatiallyFilteredImage(r+1,c+1,:), levelSpatiallyFilteredImage(r+2,c+2,:));
%                 levelGrayScaleGradientImage(r,c) = 2;%%
            end
            if (levelGrayScaleGradientImage(r,c) > 67.5 && levelGrayScaleGradientImage(r,c) <= 112.5)...
                    || (levelGrayScaleGradientImage(r,c) + 180 > 67.5 && levelGrayScaleGradientImage(r,c) + 180 <= 112.5)
                % horizontal
                spatialConsistentPixels = cat(4, levelSpatiallyFilteredImage(r,c-2,:), levelSpatiallyFilteredImage(r,c-1,:),...
                    levelSpatiallyFilteredImage(r,c,:), levelSpatiallyFilteredImage(r,c+1,:), levelSpatiallyFilteredImage(r,c+2,:));
%                 levelGrayScaleGradientImage(r,c) = 3;%%
            end
            if (levelGrayScaleGradientImage(r,c) > 112.5 && levelGrayScaleGradientImage(r,c) <= 157.5)...
                    || (levelGrayScaleGradientImage(r,c) + 180 > 112.5 && levelGrayScaleGradientImage(r,c) + 180 <= 157.5)
                % second diagonal
                spatialConsistentPixels = cat(4, levelSpatiallyFilteredImage(r-2,c+2,:), levelSpatiallyFilteredImage(r-1,c+1,:),...
                    levelSpatiallyFilteredImage(r,c,:), levelSpatiallyFilteredImage(r+1,c-1,:), levelSpatiallyFilteredImage(r+2,c-2,:));
%                 levelGrayScaleGradientImage(r,c) = 4;%%
            end

            meanVal = mean(spatialConsistentPixels, 4);
            spatialConsistentPixels = sqrt(mean(spatialConsistentPixels.^2, 3));
            spatialConsistentPixels = spatialConsistentPixels(:);
            sigmat2 = mean((spatialConsistentPixels - sqrt(mean(meanVal.^2, 3))).^2);
            sigmac2 = max(0, sigmat2 - sigma2);
            levelSpatiallyFilteredImage(r,c,:) = meanVal + sigmac2 / (sigmac2 + sigma2) * (levelSpatiallyFilteredImage(r,c,:) - meanVal);
            
            if isnan(levelSpatiallyFilteredImage(r,c,3))
                break;
            end
        end
    end
    % spatial fusion
    levelDifferenceImage = imresize(differenceImage, [size(levelSpatiallyFilteredImage, 1), size(levelSpatiallyFilteredImage, 2)]);
    refPyramid{level} = levelDifferenceImage .* levelSpatiallyFilteredImage...
        + (1 - levelDifferenceImage) .* imresize(refPyramid{level - 1}, [size(refPyramid{level},1), size(refPyramid{level},2)], 'bilinear');
    % multi-scale fusion
    % reuse the consistent pixel map for all levels
    levelConsistentPixelMapR = imresize(baseConsistentPixelMapR, [rows, cols], 'near');
    % get the set of all frames at this level
    levelImageSet = [];
    for i = 1 : length(imageSet)
        if level < length(refPyramid)
            ithImage = imresize(imageSet{i},[rows,cols], 'bilinear');
        else
            ithImage = imageSet{i};
        end
        ithImage = backwardTransform(ithImage, homographyFlowPyramidSet{level}{i});%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%
        levelImageSet = cat(4, levelImageSet, ithImage);
    end
    % use consistent image to compute mean value and variance
    levelGrayScaleImageSet = zeros(size(levelImageSet,1), size(levelImageSet,2), size(levelImageSet,4));
    for i = 1 : length(imageSet)
        levelGrayScaleImageSet(:,:,i) = rgb2gray(levelImageSet(:,:,:,i));
    end
    % get consistent image set
    levelConsistentImageSet = levelGrayScaleImageSet; %getConsistentImageSet(levelGrayScaleImageSet, homographyFlowPyramidSet{level});
    levelConsistentPixelMap = levelConsistentPixelMapR > 0;
    levelConsistentImageSet = levelConsistentImageSet .* levelConsistentPixelMap;
    omegaMap = abs(levelConsistentImageSet - repmat(rgb2gray(refPyramid{level}), [1, 1, size(levelConsistentImageSet,3)]))...
        < repmat(3 * sqrt(sigmat2MapSet{level}), [1, 1, size(levelConsistentImageSet,3)]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omegaMap = sum(omegaMap .* levelConsistentPixelMap, 3);
    omegaMap = repmat(sqrt(omegaMap / length(imageSet)), [1 1 3]);
    refPyramid{level} = omegaMap .* refPyramid{level}...
        + (1 - omegaMap) .* imresize(refPyramid{level - 1}, [size(refPyramid{level},1), size(refPyramid{level},2)], 'bilinear');
end

imwrite(uint8(refPyramid{length(refPyramid)}), [result_path,'_fixbug4.png']);
