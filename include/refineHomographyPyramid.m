function homographyPyramid = refineHomographyPyramid(homographyPyramid)
[~, layerNum] = size(homographyPyramid);
lambda = 0.1;
for level = 2:layerNum
    homographyLevel = homographyPyramid{level};
    [rows, cols] = size(homographyLevel);
    num = rows * cols;
    
    R = []; %zeros(num * 9, 1);
    D = zeros(num * 4 * 9, num * 9);
    for r = 1 : rows
        for c = 1 : cols
            homography = homographyLevel(r,c).homographies.T;
            %pointNumber = homographyLevel(r,c).pointNumber;
            %if pointNumber < 8 && level > 1
                %upperLevel = homographyPyramid{level - 1};
                %homography = upperLevel(ceil(r / 2), ceil(c / 2)).homographies.T;
            %end
            R = [R; homography(:)];
            cur_ind = (r - 1) * cols + c;
            % determine neighbors
            cur_r = r;
            cur_c = c;
            for neighbor = 1 : 4
                count = 0;
                switch neighbor
                    case 1
                        cur_r = r - 1; cur_c = c;
                    case 2
                        cur_r = r + 1; cur_c = c;
                    case 3
                        cur_r = r; cur_c = c - 1;
                    case 4
                        cur_r = r; cur_c = c + 1;
                end
                if cur_r <= rows && cur_r >= 1 && cur_c <= cols && cur_c >= 1
                    ind = (cur_r - 1) * cols + cur_c;
                    row_to_fill = (cur_ind - 1) * 4 * 9 + (neighbor - 1) * 9;
                    col_to_fill = (ind - 1) * 9;
                    for shift = 1 : 9
                        D(row_to_fill + shift, col_to_fill + shift) = -1;
                    end
                    count = 1;
                end
                if count == 1
                    row_to_fill = (cur_ind - 1) * 4 * 9 + (neighbor - 1) * 9;
                    col_to_fill = (cur_ind - 1) * 9;
                    for shift = 1 : 9
                        D(row_to_fill + shift, col_to_fill + shift) = count;
                    end
                end
            end % loop of neighborsrows
        end % loop of cols
    end % loop of rows
    
    A = eye(num * 9, num * 9) + lambda / 2 * (D' * D);
    b = R;
    
    H = A \ b;
    for r = 1 : rows
        for c = 1 : cols
            ind = (r - 1) * cols + c;
            homographyLevel(r,c).homographies.T = reshape(H((ind - 1) * 9 + 1: ind * 9), 3, 3);
        end
    end
    
    homographyPyramid{level} = homographyLevel;
end