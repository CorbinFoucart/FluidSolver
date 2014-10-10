% Triangle Returner
% Corbin Foucart
% -------------------------------------------------------------------------
% takes in a list of vertex nubmers and a triangulation connectivity list
% and returns the triangle numbers of every triangle contained in the list
% of points.

function rnTriNums = adif2d_tri_returner(region_points, conn_list)

% If conn_list contains a vertex number, keep track of it
sum_mat = zeros(size(conn_list));
for i = 1:length(region_points)  
    sum_mat = sum_mat + (conn_list == region_points(i));
end

% return indices in the matrix that sum to 3, because we want all the
% triangles that contain exactly 3 points in our point list.
ind = find((sum(sum_mat,2) == 3));
       
rnTriNums = ind;

end
