function homographyLevel = refine_core(homographyLevel, lambda)

MAX_ITER = 50;

[rows, cols] = size(homographyLevel);
num = rows * cols;

R = zeros(3, 3, num);
D = zeros(num, num);
count_record = zeros(num, 1);

% pre-define R, D and count_record, initialize H
for r = 1 : rows
    for c = 1 : cols
        cur_ind = (r-1) * cols + c;
        homography = homographyLevel(r,c).homographies.T;
        R(:,:,cur_ind) = homography;
        
        count = 0;
        for neighbor = 1 : 4
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
                D(cur_ind, ind) = 1;
                count = count + 1;
            end
        end
        count_record(cur_ind) = count;
    end
end
H = R;
obj_old = cal_object(H,R,D,lambda);
disp(obj_old);

for iter = 1 : MAX_ITER
    % update one coordinate (one H node) each time
    for i = 1 : num
        temp = sum(H(:,:,find(D(i,:))), 3);
        H(:,:,i) = (R(:,:,i)+lambda*temp) / (1+lambda*count_record(i));
    end
    
    obj_val = cal_object(H, R, D, lambda);
    if abs(obj_val-obj_old) < 0.0001
        break
    end
    obj_old = obj_val;
    disp([iter, obj_val]);
end

for r = 1 : rows
    for c = 1 : cols
        ind = (r - 1) * cols + c;
        homographyLevel(r,c).homographies.T = H(:,:,ind);
    end
end

end

function obj_val = cal_object(H, R, D, lambda)

obj_val = 0;
num = size(H, 3);
for i = 1 : num
    temp = 0;
    neighbor = find(D(i,:));
    for j = 1 : length(neighbor)
        temp = temp + norm(H(:,:,i)-H(:,:,neighbor(j)), 'fro')^2;
    end
    val = norm(H(:,:,i)-R(:,:,i), 'fro')^2 + lambda*temp;
    obj_val = obj_val + val;
end

end