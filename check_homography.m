clear all

addpath(genpath('./include'));

imageNum = 10;
ref = 5;    % set the reference image to be the 5th one

base_dir = '/localdisk/xyang/PS_data/burstimages_v1/';
name1 ='Bookshelf_2';
image_path = [base_dir, name1];


% load images
imageSet = cell(1, imageNum);
ratio = 1;
for i = 1 : imageNum
    image_dir = fullfile(image_path, [num2str(i - 1), '.jpg']);
    image_temp = rgb2gray(imresize(imread(image_dir), ratio));
    imageSet{i} = getPyramids(image_temp);
end

layerNum = length(imageSet{ref});
[rows, cols] = size(imageSet{ref}{end});

% load homography flow
version = 'fix12';
%version = 'ori';
name2 = [name1, '_hnew_', version, '.mat'];
homography_path = ['/localdisk/xyang/PS_data/', name2];
load(homography_path);

% inverse transform
TransferSet = cell(1, imageNum);
for i = 1 : length(imageSet)
    disp(i);
    if i == ref
        TransferSet{i} = imageSet{i};
        continue
    end
    
    Transfer = cell(1, layerNum);
    for j = 1 : layerNum
        ithImage = backwardTransform(imageSet{i}{j}, homographyFlowPyramidSet{j}{i});%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%
        Transfer{j} = ithImage;
    end
    TransferSet{i} = Transfer;
end

% calculate errors
error_record = zeros(layerNum, imageNum);
for i = 1 : imageNum
    for j = 1 : layerNum
        error = norm(double(TransferSet{i}{j})-double(imageSet{ref}{j}), 'fro');
        error_record(j, i) = error;
    end
end

avg_error = mean([error_record(:, 1:ref-1),error_record(:, ref+1:end)], 2)
