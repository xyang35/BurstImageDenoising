function ConsistentPixelMap = Combine_strategy(RefConsistentPixelMap, MedConsistentPixelMap, ref)
% Combinign the consistent pixel selection result from reference based and median based method

[rows, cols, imageNum] = size(RefConsistentPixelMap);
ConsistentPixelMap = zeros(size(RefConsistentPixelMap));

undecided = ones(rows, cols);

% case 1: union of the two
for r = 1 : rows
    for c = 1 : cols
        if MedConsistentPixelMap(r,c,ref) == 1
            ConsistentPixelMap(r,c,:) = RefConsistentPixelMap(r,c,:) | MedConsistentPixelMap(r,c,:);
            undecided(r,c) = 0;
        end
    end
end

reliableNumber = floor(imageNum / 2);
consistentPixelNumMap = sum(MedConsistentPixelMap, 3);
consistentPixelNumMap = consistentPixelNumMap > reliableNumber;

% find connected component of undecided pixels
CC = bwconncomp(undecided);

% majority voting for each connected omponent
for i = 1 : CC.NumObjects
    temp = consistentPixelNumMap(CC.PixelIdxList{i});
    majority = mode(double(temp));
    consistentPixelNumMap(CC.PixelIdxList{i}) = majority;
end

for r = 1 : rows
    for c = 1 : cols
        if undecided(r,c) == 1
            if consistentPixelNumMap(r,c) == 1
                % case 2: median based result is reliable
                ConsistentPixelMap(r,c,:) = MedConsistentPixelMap(r,c,:);
            else
                % case 3: median based result is not reliable
                ConsistentPixelMap(r,c,:) = RefConsistentPixelMap(r,c,:);
            end
        end
    end
end

% morphological filter
for i = 1 : imageNum
    ConsistentPixelMap(:, :, i) = bwmorph(ConsistentPixelMap(:,:,i), 'majority');
end
