% Rectangular Random Voronoi Mesh Generator
% Corbin Foucart
% ------------------------------------------------------------------------
% For a rectangular region [a,b]x[c,d]

clear all; close all; clc;

% ------------------------------- Input ------------------------------- %
% input data for each dimension
a = 0;
b = 1;

c = 0;
d = 1;

N_boundary_points = 20;
N_inner_points = N_boundary_points^2;
N_corner_points = 4;

% For saving to a data file later
filename = 'zVrm20';

% Generate vertex points on the domain with our point generation function
% sf: scaling factor between 0 and 0.5, doesn't place points that
% percentage away from the interval wall. 

sf = 0;
points = vg2d_point_creation(a,b,c,d,N_boundary_points,...
    N_inner_points, N_corner_points, sf);


% ----------------------- Delaunay Triangulation ------------------------%

tri = delaunay(points);
dtri = delaunayTriangulation(points);
dtri_points = dtri.Points;              % row is vertex num
conn_list = dtri.ConnectivityList;      % row is triangle num
CC = circumcenter(dtri);                % row is triangle num

% ------------------------- Property Tracking  --------------------------%

% Create sparse neighbor matrix to relate each vertex number to the
% vertices of its neighbors. 
% rows: vertex number, columns: neighbor vertex numbers, 1 if true
neighbor_mat = vg2d_neigh_mat_creation(dtri_points, conn_list);

% Create edge matrix 
% rows: edge number (arbitrarily), columns: 2 vertex numbers of the edge
edge_matrix = vg2d_edge_mat_creation(neighbor_mat);

% Create matrix giving two adjacent triangles from an edge
% row is the edge num, columns are triangle numbers
adj_triangles = vg2d_adj_triangles_list_creation(edge_matrix, conn_list);
 
% ------------------------- Voronoi Diagram -----------------------------%
% Point here is to connect lines

% vertex_vlines: sparse matrix to relate vertex number (row) to voronoi line (column)
% v_lines: row is voronoi line columns are points, 5th column will be length.
% [x1, y1, x2, y2]

[v_lines, vertex_vlines] = vg2d_voronoi_lines_creation(edge_matrix,...
    dtri_points, CC, adj_triangles);

% ------------------------ Voronoi Geometric Data --------------------- %

% Mark which points are boundary points, append to dtri_points
dtri_points = vg2d_boundary_checker(a,b,c,d, dtri_points);

% Get lengths of voronoi segments, append to v_lines
vdistances = zeros(length(v_lines(:,1)),1);
for i = 1:length(v_lines(:,1))
    vdistances(i) = sqrt((v_lines(i,3)-v_lines(i,1))^2 ...
        + (v_lines(i,4)-v_lines(i,2))^2);
end
% 5th column distance of the voronoi lines
v_lines = [v_lines, vdistances]; 

% important to note that edge number is also voronoi line number. 
% 3rd column of edge matrix contains length of edge
edge_matrix = [edge_matrix, vdistances];


% Areas, barycenters of voronoi boxes
% Row is voronoi vertx num, column is area/ [x,y] of barycenter.
 [voronoi_areas, vertex_barycenters] = ...
    vg2d_area_barycenter_calc(vertex_vlines, v_lines, dtri_points);

% ---------------------- Plotting, Visualization ----------------------- %
figure()
hold all
scatter(points(:,1), points(:,2), 'ko')
scatter(CC(:,1), CC(:,2), 'k.')
% get handle to scattergroup object
h = get(gca,'children');
% change size of markers 
set(h, 'sizedata', 50)
for i = 1:length(v_lines)
    plot(line([v_lines(i,1), v_lines(i,3)], [v_lines(i,2), v_lines(i,4)], ... 
       'color', 'r', 'LineWidth', 4))
end
axis([a b c d])
axis square
triplot(dtri, 'color', [0.5 0.5 0.5], 'LineWidth', 0.5, 'LineStyle', '--')

figure()
voronoi(dtri)

% Save Meshing to file to be used by solving programs later
 save(filename);


% Gravevyard Code
%{
while(true)
    % prob_cc: rows are triangle_num, columns are coordinates
    prob_cc = zeros(length(CC(:,1)),2);
    repair_pts = zeros(length(CC(:,1)),2); % must have 0s removed
    for i = 1:length(CC(:,1))
        cct = CC(i,:);
        % check if points are outside boundaries, create new midpoints
        if cct(1) < a 
            prob_cc(i,:) = CC(i,:);
            prob_pts =  sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 2);
            repair_pts(i,:) = [mean(prob_pts(1:2, 2)), mean(prob_pts(1:2, 3))];
            prob_pt_ind = prob_pts(1:2,1)'; % select indices of points closest to boundary

        end
        if cct(1) > b
            prob_cc(i,:) = CC(i,:);
            prob_pts =  flipud(sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 2));
            repair_pts(i,:) = [mean(prob_pts(1:2, 2)), mean(prob_pts(1:2, 3))];
        end
        if cct(2) < c
            prob_cc(i,:) = CC(i,:);
            prob_pts =  sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 1);
            repair_pts(i,:) = [mean(prob_pts(1:2, 2)), mean(prob_pts(1:2, 3))];
        end
        if cct(2) > d
            prob_cc(i,:) = CC(i,:);
            prob_pts =  flipud(sortrows([conn_list(i,:)', dtri_points(conn_list(i,:),:)], 1));
            repair_pts(i,:) = [mean(prob_pts(1:2, 2)), mean(prob_pts(1:2, 3))];
        end
    end

    if prob_cc == zeros(length(CC(:,1)),2)
        break
    end

    % add repair points, remove redundancies
    points = unique([points; repair_pts],'rows');

    % re triangulate 
    tri = delaunay(points);
    dtri = delaunayTriangulation(points);
    dtri_points = dtri.Points;              % row is vertex num
    conn_list = dtri.ConnectivityList;      % row is triangle num
    CC = circumcenter(dtri);                % row is triangle num
    [mcon ncon] = size(conn_list);
end

%}
