% Voronoi 2D mesh generator, Voronoi Areas, Barycenters
% Corbin Foucart
%-------------------------------------------------------------------------
% Calculates, areas, barycenters of each voronoi box

function [rn_va, rn_vb] = vg2d_area_barycenter_calc(...
    vertex_vlines, v_lines, dtri_points)

% Areas of voronoi boxes
% Row is voronoi vertx num, column is area.
voronoi_areas = zeros(length(vertex_vlines(:,1)),1);

% Vertex barycenters, row is vertex nubmer, columns 1,2 are x and y 
% coordinates of the barycenter corresponding to the vertex number. 
vertex_barycenters = zeros(length(vertex_vlines(:,1)),1);

% Computation
for i = 1:length(vertex_vlines(:,1))
    ind = find(vertex_vlines(i,:));
    poly_pts = [];
    for j = 1:length(ind)
        poly_pts = [poly_pts; v_lines(ind(j),1:2); v_lines(ind(j), 3:4)];
    end
    if dtri_points(i,3) == 1 % include vertex point in polygon if boundary
        poly_pts = [poly_pts; dtri_points(i,1:2)];
    end
    poly_pts = unique(poly_pts,'rows');
    voronoi_areas(i) = p_area_calc(poly_pts(:,1), poly_pts(:,2));
    vertex_barycenters(i,1:4) = barycenter_calc(poly_pts(:,1), poly_pts(:,2));
end

vertex_barycenters = vertex_barycenters(:,2:3);

rn_va = voronoi_areas;
rn_vb = vertex_barycenters;
end

