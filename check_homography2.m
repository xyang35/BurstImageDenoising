clear all

addpath(genpath('./include'));

imageNum = 10;
ref = 5;    % set the reference image to be the 5th one

base_dir = '/localdisk/xyang/PS_data/burstimages_v1/';
name1 = 'Bookshelf_2';
image_path = [base_dir, name1];


% load images
imageSet = cell(1, imageNum);
ratio = 1;
for i = 1 : imageNum
    image_dir = fullfile(image_path, [num2str(i - 1), '.jpg']);
    imageSet{i} = rgb2gray(imresize(imread(image_dir), ratio));
end

[rows, cols] = size(imageSet{5});

% load homography flow
version = 'fix12';
%version = 'ori';
name2 = [name1, '_hnew_', version, '.mat'];
homography_path = ['/localdisk/xyang/PS_data/', name2];
load(homography_path);

% inverse transform
TransferSet = uint8(zeros(rows, cols, imageNum));
for i = 1 : length(imageSet)
    disp(i);
    if i == ref
        TransferSet(:,:,i) = imageSet{i};
    end
    ithImage = backwardTransform(imageSet{i}, homographyflow{i});%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%ATTENTION HERE%%%%%
    TransferSet(:,:,i) = ithImage;
end

% calculate errors
error_record = zeros(1, imageNum);
for i = 1 : imageNum
    error = norm(double(TransferSet(:,:,i))-double(imageSet{ref}), 'fro');
    error_record(i) = error;
end

avg_error = mean([error_record(1:ref-1),error_record(ref+1:end)])
