% Voronoi 2D mesh generator, edge matrix generator
% Corbin Foucart
%-------------------------------------------------------------------------
% Takes in the neighbor matrix and creates an 'edge matrix'
% This is a matrix where the rows number the edges, and the 2 columns and
% contain the two vertex numbers that form the edge. 

% row: edge number, columns: two vertex numbers

function return_edge_mat = vg2d_edge_mat_creation(neighbor_mat)

% There are redunancies in edge information in the neighbor matrix,
% so take once an edge is recorded, remove the corresponding redundancy
% from the matrix that will occur later. Remove empty entries.
bkkping_neighbor_mat = neighbor_mat;
edge_matrix = zeros(sum(neighbor_mat(:)), 2);
counter = 1;
for i = 1:length(neighbor_mat(:,1))
    for j = 1:length(neighbor_mat(:,1))
        if bkkping_neighbor_mat(i,j) ~= 0
            bkkping_neighbor_mat(j,i) = 0; % removing the redundancy
            edge_matrix(counter,:) = [i,j];
            counter = counter + 1;
        end
    end
end

% remove empty entries
return_edge_mat = edge_matrix(any(edge_matrix,2),:); 
end