% Triangle Returner
% Corbin Foucart
% -------------------------------------------------------------------------
% takes in a list of vertex numbers and a edge matrix list
% and returns the edge numbers of every edge contained in the list
% of points.

function rnEdgeNums = adif2d_edge_returner(region_points, edge_matrix)

edge_matrix = edge_matrix(:,1:2);
% If conn_list contains a vertex number, keep track of it
sum_mat = zeros(size(edge_matrix));
for i = 1:length(region_points)  
    sum_mat = sum_mat + (edge_matrix == region_points(i));
end

% return indices in the matrix that sum to 2, because we want all the
% edges that contain exactly 2 points in our point list.
ind = find((sum(sum_mat,2) == 2));
       
rnEdgeNums = ind;
end
