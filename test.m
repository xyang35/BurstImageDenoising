% %close all; clear all;
% addpath(genpath('./include'));
% 
% %% Initialization ad global parameters
% 
% imageNum = 10;
% ref = 5;    % set the reference image to be the 5th one
% 
% %base_dir = '/localdisk/xyang/PS_data/burstimages_v1/';
% base_dir = './';
% name = 'Building_1';
% path = [base_dir, name];
% 
% %result_path = ['/localdisk/xyang/PS_data/', name];
% result_path = ['result/', name];
% method = 'nofix1';
% 
% 
% % import burst of images
% imageSet = cell(1, imageNum);
% ratio = 1;
% for i = 1 : imageNum
%     image_dir = fullfile(path, [num2str(i - 1), '.jpg']);
%     imageSet{i} = imresize(imread(image_dir), ratio);
%     imageSet{i} = imageSet{i}(400:1300, :, :);
%     imwrite(imageSet{i}, fullfile(path, [num2str(i-1), '_crop.jpg']), 'jpg');
% end
% refImage = imageSet{ref}; % get reference image
% refPyramid = getPyramids(refImage); % get pyramid for reference image
% 
% % scale level (s=1) is obtained after one down-sampling operation
% PyHeight = length(refPyramid);
% FEATURELEVEL = PyHeight - 1;
% CONSISTLEVEL = PyHeight - 1;
% [refFeatures, refPoints] = getFeatures(refPyramid{FEATURELEVEL}, FEATURELEVEL); % get features and valid points of reference image
% 
% %% Homography flows
% 
% % homography flows for all base images (base: coarsest level of the pyramid)
% homographyFlowPyramidSet = cell(1, length(refPyramid));
% ImagePyramidSet = cell(1, imageNum);
% ImagePyramidSet{ref} = refPyramid;
% for j = 1 : length(refPyramid)
%     homographyFlowPyramidSet{j} = cell(1,length(imageSet));
% end
% % record scale 1 gray scale image set for consistent pixel selection
% Scale1GrayImageSet = zeros(size(refPyramid{CONSISTLEVEL}, 1), size(refPyramid{CONSISTLEVEL}, 2), length(imageSet));
% 
% for i = 1 : length(imageSet)
%      if i == ref
%          Scale1GrayImageSet(:,:,i) = rgb2gray(refPyramid{CONSISTLEVEL});
%          continue;
%      end
%      pyramid = getPyramids(imageSet{i});
%      [homographyFlowPyramid, homographyPyramid] = getHomographyFlowPyramidWithRefFeatures(refPyramid, refFeatures, refPoints, pyramid, FEATURELEVEL);
%      for j = 1 : length(homographyFlowPyramid)
%          homographyFlowPyramidSet{j}{i} = homographyFlowPyramid{j};
%      end
%      Scale1GrayImageSet(:,:,i) = rgb2gray(pyramid{CONSISTLEVEL});
%      % just for debugging record
%      %ImagePyramidSet{i} = pyramid;
%      disp(['homography ', num2str(i), ' complete']);
% end
% 
% %load([result_path, '_', method, '.mat']);
% %homographyflow = homographyFlowPyramidSet{end};
% save(['result/', name, '_',method,'.mat'], 'homographyFlowPyramidSet', '-v7.3');
% %return
% 
% 
% %% Consistent pixel
% tic
% disp('Consisten pixel selection ...');
% tau = 10; % threshold for selecting consistent pixels
% 
% % compute median image of all consistent level images (for median value based consistent pixel selection)
% ConsistentImageSet = getConsistentImageSet(Scale1GrayImageSet, homographyFlowPyramidSet{CONSISTLEVEL});
% MedianImage = median(ConsistentImageSet, 3);
% IntegralMedImage = integralImage(MedianImage); % integral image of the median image
% % integral images of all consistent level images
% IntegralImageSet = zeros(size(refPyramid{CONSISTLEVEL}, 1) + 1, size(refPyramid{CONSISTLEVEL}, 2) + 1, length(imageSet));
% for i = 1 : length(imageSet)
%     IntegralImageSet(:, :, i) = integralImage(Scale1GrayImageSet(:, :, i));
% end
% 
% % integral images of all consistent images
% ConsistentIntegralImageSet = zeros(size(ConsistentImageSet, 1) + 1, size(ConsistentImageSet, 2) + 1, size(ConsistentImageSet, 3));
% for i = 1 : size(ConsistentIntegralImageSet,3)
%     ConsistentIntegralImageSet(:, :, i) = integralImage(ConsistentImageSet(:, :, i));
% end
% 
% % set of consistent pixel indexes: reference based and median based
% RefConsistentPixelMap = zeros(size(Scale1GrayImageSet));
% MedConsistentPixelMap = zeros(size(Scale1GrayImageSet));
% halfWidth = 2; halfHeight = 2;    % for patch implementation
% rows = size(RefConsistentPixelMap, 1);
% cols = size(RefConsistentPixelMap, 2);
% 
% % rs = 1 : 10 : rows; cs = 1 : 10 : cols;
% % [x, y] = meshgrid(cs, rs);
% % quiver(x,y,baseHomographySet(rs,cs,1,i),baseHomographySet(rs,cs,2,i));
% 
% 
% for r = 1 : rows
%     for c = 1 : cols
%         % first get the pixel from the reference image
%         sR = max(1, r - halfHeight);
%         sC = max(1, c - halfWidth);
%         eR = min(rows, r + halfHeight);
%         eC = min(cols, c + halfWidth);
%         pixNum = (eR - sR + 1) * (eC - sC + 1);
%         refPix = IntegralImageSet(eR+1,eC+1,ref) - IntegralImageSet(eR+1,sC,ref) - IntegralImageSet(sR,eC+1,ref) + IntegralImageSet(sR,sC,ref);
%         refPix = refPix / pixNum;
%         medPix = IntegralMedImage(eR+1,eC+1) - IntegralMedImage(eR+1,sC) - IntegralMedImage(sR,eC+1) + IntegralMedImage(sR,sC);
%         medPix = medPix / pixNum;
% 
%         % reference-based
%         RefConsistentPixelMap(r,c,ref) = 1;
%         % from reference to left
%         for i = ref - 1 : -1 : 1
%             % make use of ConsistenImage to get consistent pixels
%             iPix = ConsistentIntegralImageSet(eR+1,eC+1,i) - ConsistentIntegralImageSet(eR+1,sC,i) - ConsistentIntegralImageSet(sR,eC+1,i) + ConsistentIntegralImageSet(sR,sC,i);
%             iPix = iPix / pixNum;
%             if abs(refPix - iPix) < tau
%                 RefConsistentPixelMap(r,c,i) = 1;
%             else
%                 break;
%             end
%         end
%         % from reference to right
%         for i = ref + 1 : length(imageSet)
%             iPix = ConsistentIntegralImageSet(eR+1,eC+1,i) - ConsistentIntegralImageSet(eR+1,sC,i) - ConsistentIntegralImageSet(sR,eC+1,i) + ConsistentIntegralImageSet(sR,sC,i);
%             iPix = iPix / pixNum;
%             if abs(refPix - iPix) < tau
%                 RefConsistentPixelMap(r,c,i) = 1;
%             else
%                 break;
%             end
%         end
% 
%         % median-based
%         for i = 1 : length(imageSet)
%             iPix = ConsistentIntegralImageSet(eR+1,eC+1,i) - ConsistentIntegralImageSet(eR+1,sC,i) - ConsistentIntegralImageSet(sR,eC+1,i) + ConsistentIntegralImageSet(sR,sC,i);
%             iPix = iPix / pixNum;
%             if abs(medPix - iPix) < tau
%                 MedConsistentPixelMap(r,c,i) = 1;
%             end
%         end
%     end
% end
% 
% % combine median consistent pixels and reference consistent pixels
% 
% disp('Combine strategy ...');
% ConsistentPixelMap = Combine_strategy(RefConsistentPixelMap, MedConsistentPixelMap, ref);
% %ConsistentPixelMap = zeros(size(RefConsistentPixelMap));
% %reliableNumber = floor(imageNum / 2);
% %consistentPixelNumMap = sum(MedConsistentPixelMap, 3);
% %consistentPixelNumMap = consistentPixelNumMap > reliableNumber;
% %% perform majority filter
% %consistentPixelNumMap = bwmorph(consistentPixelNumMap, 'majority');
% %
% %for r = 1 : rows
% %    for c = 1 : cols
% %        % case 1: union of the two
% %        if MedConsistentPixelMap(r,c,ref) == 1
% %            ConsistentPixelMap(r,c,:) = RefConsistentPixelMap(r,c,:) | MedConsistentPixelMap(r,c,:);
% %        % case 2: judge if median based is reliable
% %        elseif consistentPixelNumMap(r,c) == 1
% %            % median based result is reliable
% %            ConsistentPixelMap(r,c,:) = MedConsistentPixelMap(r,c,:);
% %        else
% %            % median based result is not reliable
% %            ConsistentPixelMap(r,c,:) = RefConsistentPixelMap(r,c,:);
% %        end
% %    end
% %end
% 
% sumConsistentPixelMap = sum(ConsistentPixelMap, 3);
% [rs, cs] = find(sumConsistentPixelMap == 0);
% ConsistentPixelMap(rs, cs, ref) = 1;
% 
% % reuse the consistent pixel map for all levels by upsampling or downsampling
% AllConsistentPixelMap = cell(1, length(refPyramid));
% for level = 1 : length(refPyramid)
%     rows = size(refPyramid{level}, 1);
%     cols = size(refPyramid{level}, 2);
%     if level == CONSISTLEVEL
%         AllConsistentPixelMap{level} = ConsistentPixelMap;
%     else
%         AllConsistentPixelMap{level} = imresize(ConsistentPixelMap, [rows, cols], 'near');
%     end
% end
% 
% %% fusion stage
% 
% disp('Fusion stage')
% % first estimate noise
% fineRefImage = refPyramid{length(refPyramid)};
% fineGrayScaleRefImage = double(rgb2gray(fineRefImage));
% nontextureMap = imresize(edge(MedianImage), size(fineGrayScaleRefImage), 'near');
% inds = find(nontextureMap == 0);
% fineMedianImage = imresize(MedianImage, size(fineGrayScaleRefImage), 'bilinear');
% sigma2 = computeSigma2FromDiffVector(fineGrayScaleRefImage(inds) - fineMedianImage(inds));
% 
% % second perform temporal fusion
% disp('temporal fusion...')
% sigmat2MapSet = cell(1,length(refPyramid));
% levelImageSet_record = cell(1, length(refPyramid));
% for level = 1 : length(refPyramid)
%     rows = size(refPyramid{level}, 1);
%     cols = size(refPyramid{level}, 2);
%     levelConsistentPixelMap = AllConsistentPixelMap{level};
% 
%     % get the set of all frames at this level
%     levelImageSet = [];
%     for i = 1 : length(imageSet)
%         if i == ref
%             levelImageSet = cat(4, levelImageSet, refPyramid{level});
%             continue;
%         else
%             if level < length(refPyramid)
%                 ithImage = imresize(imageSet{i},[rows,cols], 'bilinear');
%             else
%                 ithImage = imageSet{i};
%             end
%             ithImage = backwardTransform(ithImage, homographyFlowPyramidSet{level}{i});%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%
%             levelImageSet = cat(4, levelImageSet, ithImage);
%         end
%     end
%     levelImageSet_record{level} = levelImageSet;
% 
%     % use consistent image to compute mean value and variance
%     levelConsistentImageSet = zeros(size(levelImageSet,1), size(levelImageSet,2), size(levelImageSet,4));
%     for i = 1 : length(imageSet)
%         levelConsistentImageSet(:,:,i) = rgb2gray(levelImageSet(:,:,:,i));
%     end
%     % get consistent image set
%     levelConsistentPixelMap = levelConsistentPixelMap > 0;
%     levelConsistentImageSet = levelConsistentImageSet .* levelConsistentPixelMap;
%     meanImage = sum(levelConsistentImageSet, 3) ./ sum(double(levelConsistentPixelMap), 3);
%     % sigma_t^2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     sigmat2MapSet{level} = sum((levelConsistentImageSet - repmat(meanImage, [1 1 size(levelConsistentImageSet, 3)])) .^ 2, 3) ./ sum(levelConsistentPixelMap, 3);
%     % sigma_c^2
%     sigmac2Map = sigmat2MapSet{level} - sigma2;
%     sigmac2Map = (sigmac2Map > 0) .* sigmac2Map;
%     sigmac2Map = sigmac2Map ./ (sigmac2Map + sigma2);
%     sigmac2Map = repmat(sigmac2Map, [1 1 3]);
%     levelConsistentPixelMap = reshape(levelConsistentPixelMap,...
%         [size(levelConsistentPixelMap,1), size(levelConsistentPixelMap,2), 1, size(levelConsistentPixelMap,3)]);
%     levelConsistentPixelMap = repmat(levelConsistentPixelMap, [1,1,3,1]);
%     meanImage = sum(double(levelImageSet) .* double(levelConsistentPixelMap), 4) ./ sum(double(levelConsistentPixelMap), 4);
%     refPyramid{level} = meanImage + sigmac2Map .* (double(refPyramid{level}) - meanImage);
% end
% %figure;imshow(uint8(refPyramid{length(refPyramid)}))
% imwrite(uint8(refPyramid{length(refPyramid)}), [result_path,'_temporalonly_', method, '.png']);
% % snapshot of current state
% save([result_path,'_',method,'_temporal_state.mat'])

load('result/Bookshelf_2_nofix1_temporal_state.mat')
figure;imshow(uint8(refPyramid{length(refPyramid)}))
method='fix8';
% compute textureness probability p_tex
disp('spatial fusion...')
differenceImage = [];    % absolute difference between the pixel and its 4 neighbors
rs = [1 2 2 3];
cs = [2 1 3 2];
for i = 1 : 4
    h = zeros(3,3);
    h(2,2) = 1;
    h(rs(i), cs(i)) = -1;
    differenceImage = cat(3, differenceImage, abs(conv2(fineGrayScaleRefImage, h, 'same')));
end
differenceImage = max(differenceImage, [], 3);
p_tex = 1 ./ (1 + exp(-5 * (differenceImage / sqrt(sigma2) - 3)));
figure; image((p_tex>0.01)*255);
%p_tex(:,:) = 0;    % for testing spatial fusion
p_tex = repmat(p_tex, [1 1 3]);
for level = 2 : length(refPyramid)
    levelSpatiallyFilteredImage = refPyramid{level};
    p_tex_level = imresize(p_tex, [size(levelSpatiallyFilteredImage, 1), size(levelSpatiallyFilteredImage, 2)]);
    
    [rows, cols, ~] = size(levelSpatiallyFilteredImage);
    formerlayerImage = imresize(refPyramid{level - 1}, [size(refPyramid{level},1), size(refPyramid{level},2)], 'bilinear');
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
            if p_tex_level(r, c, 1) > 0.01
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
                    disp(['NaN error!']);
                    return
                end
                
                % fusion
                refPyramid{level}(r,c,:) = p_tex_level(r,c,:) .* levelSpatiallyFilteredImage(r,c,:) + ...
                    (1 - p_tex_level(r,c,:)) .* formerlayerImage(r,c,:);
            end
        end
    end

    %levelDifferenceImage = imresize(p_tex, [size(levelSpatiallyFilteredImage, 1), size(levelSpatiallyFilteredImage, 2)]);
    %refPyramid{level} = levelDifferenceImage .* levelSpatiallyFilteredImage...
    %    + (1 - levelDifferenceImage) .* imresize(refPyramid{level - 1}, [size(refPyramid{level},1), size(refPyramid{level},2)], 'bilinear');

    % multi-scale fusion
    % reuse the consistent pixel map for all levels
    levelConsistentPixelMap = AllConsistentPixelMap{level};
    % get the set of all frames at this level
    levelImageSet = levelImageSet_record{level};

    % use consistent image to compute mean value and variance
    levelConsistentImageSet = zeros(size(levelImageSet,1), size(levelImageSet,2), size(levelImageSet,4));
    for i = 1 : length(imageSet)
        levelConsistentImageSet(:,:,i) = rgb2gray(levelImageSet(:,:,:,i));
    end
    % get consistent image set
    levelConsistentPixelMap = levelConsistentPixelMap > 0;
    levelConsistentImageSet = levelConsistentImageSet .* levelConsistentPixelMap;
    omegaMap = abs(levelConsistentImageSet - repmat(double(rgb2gray(uint8(refPyramid{level}))), [1, 1, size(levelConsistentImageSet,3)]))...
        < repmat(3 * sqrt(sigmat2MapSet{level}), [1, 1, size(levelConsistentImageSet,3)]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    omegaMap = sum(omegaMap .* levelConsistentPixelMap, 3);
    omegaMap = repmat(sqrt(omegaMap / length(imageSet)), [1 1 3]);
    
    %omegaMap(:,:,:) = 0.6;    % for testing multi-scale fusion
    refPyramid{level} = omegaMap .* refPyramid{level}...
        + (1 - omegaMap) .* formerlayerImage;
end
figure;imshow(uint8(refPyramid{length(refPyramid)}))
%toc
imwrite(uint8(refPyramid{length(refPyramid)}), [result_path,'_', method, '.png']);
