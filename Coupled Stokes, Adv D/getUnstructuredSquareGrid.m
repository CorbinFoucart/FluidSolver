function [nx, ny] = getUnstructuredSquareGrid(n)
% getUnstructuredSquareGrid on square [0, 1]^2

%nx = [0  0  1/2 3/10 35/100 7/10 65/100 1  1  0 1/2 1]';
%ny = [0 1/2  0  3/10 65/100 7/10 35/100 0 1/2 1  1  1]';

nv = (2 * n + 1) * (n + 1) + (n + 1) * n + 4 * n^2;

nx = zeros(nv, 1);
ny = zeros(nv, 1);

cnt = 1;
% define vertices on vertical edges
for i=1:n+1
    for j=1:2 * n + 1
        nx(cnt, 1) = (i-1) / n;
        ny(cnt, 1) = (j-1) / (2 * n);
        
        cnt = cnt + 1;
    end
end

% define vertices on horizontal edges
for i=1:n
    for j=1:n+1
        nx(cnt, 1) = 1 / (2 * n) + (i - 1) / n;
        ny(cnt, 1) = (j - 1) / n;
        
        cnt = cnt + 1;
    end
end

% define vertices in the interior of square parts
for i=1:n
    for j=1:n
        nx(cnt, 1) = (i-1) / n + 0.3 / n;
        ny(cnt, 1) = (j-1) / n + 0.3 / n;
        
        nx(cnt + 1, 1) = (i-1) / n + 0.35 / n;
        ny(cnt + 1, 1) = (j-1) / n + 0.65 / n;
        
        nx(cnt + 2, 1) = (i-1) / n + 0.7 / n;
        ny(cnt + 2, 1) = (j-1) / n + 0.7 / n;

        nx(cnt + 3, 1) = (i-1) / n + 0.65 / n;
        ny(cnt + 3, 1) = (j-1) / n + 0.35 / n;
        
        cnt = cnt + 4;
    end
end

end

