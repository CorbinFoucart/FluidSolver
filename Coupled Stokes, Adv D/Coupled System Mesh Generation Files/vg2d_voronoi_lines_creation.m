% Voronoi 2D mesh generator, voronoi line creator
% Corbin Foucart
%-------------------------------------------------------------------------
% % vertex_vlines: sparse matrix to relate vertex number (row) to voronoi line (column)
% v_lines: row is voronoi line columns are points, 5th column will be length.
% [x1, y1, x2, y2]

function [rn_v_lines, rn_vertex_vlines] = vg2d_voronoi_lines_creation(...
    edge_matrix, dtri_points, CC, adj_triangles)

% Initialize arrays for both return values

% v_lines, row is voronoi line columns are points, 5th column will be length.
v_lines = zeros(length(adj_triangles), 4);  % [x1, y1, x2, y2]
% vertex_vlines Sparse matrix to relate vertex number (row) to voronoi line (column)
vertex_vlines = zeros(length(dtri_points), length(adj_triangles));

% Calculate values
    for i = 1:length(adj_triangles)
        if adj_triangles(i,2) ~= 0
            v_lines(i,:) = [CC(adj_triangles(i,1),:), ...
                CC(adj_triangles(i,2),:)];
            bp = edge_matrix(i,:);              % get the vertex numbers
            vertex_vlines(bp(1), i) = 1;
            vertex_vlines(bp(2), i) = 1;

        end
        if adj_triangles(i,2) == 0              % boundary triangle case
            bp = edge_matrix(i,:);              % get the vertex numbers
            bv = [dtri_points(bp(1),:), dtri_points(bp(2), :)]; 
            bmp = [1/2*(bv(1)+bv(3)), 1/2*(bv(2)+bv(4))];
            v_lines(i,:) = [CC(adj_triangles(i,1),:), bmp];
            vertex_vlines(bp(1), i) = 1;
            vertex_vlines(bp(2), i) = 1;
        end
    end
    
rn_v_lines = v_lines;
rn_vertex_vlines = vertex_vlines;
    
end