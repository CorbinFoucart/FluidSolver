function [nx, ny] = getStructuredSquareGrid(n)
% getStructuredSquareGrid on square [0, 1]^2

nx = zeros((n+1)^2, 1);
ny = zeros((n+1)^2, 1);

for i=1:n+1
    for j=1:n+1
        nx((i-1) * (n+1) + j, 1) = (i-1) / n;
        ny((i-1) * (n+1) + j, 1) = (j-1) / n;
   end
end

end

