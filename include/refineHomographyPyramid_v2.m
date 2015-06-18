function homographyPyramid = refineHomographyPyramid_v2(homographyPyramid)

[~, layerNum] = size(homographyPyramid);
lambda = 0.1;

for level = 2 : layerNum
    homographyLevel = homographyPyramid{level};
    
    homographyPyramid{level} = refine_core(homographyLevel, lambda);
end 