refImage = ImagePyramidSet{1}{1};

rows = [1,size(ImagePyramidSet{1}{1},1)];
cols = [1,size(ImagePyramidSet{1}{1},2)];

homography = homographyPyramid{1}.homographies.T;
homographyFLow = zeros(rows(2) - rows(1) + 1, cols(2) - cols(1) + 1, 2);
%
x = cols(1) : cols(2);
y = rows(1) : rows(2);

for i = 1 : length(x)
    for j = 1 : length(y)
        pt = [x(i), y(j), 1];
        pt_ = pt * homography;
        pt_ = pt_ / pt_(3);
        homographyFLow(j, i, 1) = pt_(1) - pt(1);
        homographyFLow(j, i, 2) = pt_(2) - pt(2);
    end
end

refImage_transfer = forwardTransform(refImage, homographyFLow);
refImage_warp = imwarp(refImage, homographyPyramid{1}.homographies);

figure;
subplot(2,1,1);imshow(refImage_transfer);
subplot(2,1,2);imshow(refImage_warp);

refImage = ImagePyramidSet{1}{3};
homographyFLow = homographyFlowPyramidSet{3}{2};
refImage_transfer = forwardTransform(refImage, homographyFLow);
figure;imshow(refImage_transfer)

tic
tranImage = imageSet{2};
homographyFLow = homographyFlowPyramidSet{3}{2};
tranImage_transfer = backwardTransform(tranImage, homographyFLow);
figure;imshow(tranImage_transfer)
toc