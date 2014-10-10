function [nx, ny] = getRandomSquareGrid(n)
% getRandomSquareGrid on square [0, 1]^2
rng default

nNodes = (n+1)^2;

nx = zeros(nNodes, 1);
ny = zeros(nNodes, 1);

nx(1, 1) = 0;
ny(1, 1) = 0;

nx(2, 1) = 1;
ny(2, 1) = 0;

nx(3, 1) = 1;
ny(3, 1) = 1;

nx(4, 1) = 0;
ny(4, 1) = 1;

setNodes = 4;

for i=1:n-1
    setNodes = setNodes + 1;
    nx(setNodes, 1) = rand(1);
    ny(setNodes, 1) = 0;
end

for i=1:n-1
    setNodes = setNodes + 1;
    nx(setNodes, 1) = 1;
    ny(setNodes, 1) = rand(1);
end

for i=1:n-1
    setNodes = setNodes + 1;
    nx(setNodes, 1) = rand(1);
    ny(setNodes, 1) = 1;
end

for i=1:n-1
    setNodes = setNodes + 1;
    nx(setNodes, 1) = 0;
    ny(setNodes, 1) = rand(1);
end

for i=setNodes+1:nNodes
    setNodes = setNodes + 1;
    nx(setNodes, 1) = rand(1);
    ny(setNodes, 1) = rand(1);
end

end

