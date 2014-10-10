% Voronoi 2D mesh generator, neighbor matrix generator
% Corbin Foucart
%-------------------------------------------------------------------------
% Takes in triangluation points and the associated connectivity list, and
% outputs a sparse matrix detailing which vertices are neighbors.
% 
% Row of matrix: vertex number, Columns: 1's in indices for vertex numbers
% that are neighbors

function return_neigh_mat = vg2d_neigh_mat_creation(dtri_points, conn_list)

[mcon ncon] = size(conn_list);
neighbor_mat = zeros(length(dtri_points));
for i = 1:mcon
    tri_piece = conn_list(i,:);
    for j = 1:3
        point = tri_piece(j);
        for k = 1:3
            if k ~= j
                neighbor_num = tri_piece(k);
                neighbor_mat(point,neighbor_num) = 1;                
            end
        end
    end
end

return_neigh_mat = neighbor_mat;

end