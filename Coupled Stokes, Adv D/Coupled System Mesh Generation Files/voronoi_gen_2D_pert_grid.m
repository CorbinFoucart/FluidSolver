% Rectangular Voronoi Mesh Generator
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

filename = 'voronoi_pertgrid_data_10pts'
N_side_pts1 = 11;
N_side_pts2 = round((d-c)/(b-a)*N_side_pts1);

% Box Spacing
h1 = abs(b-a)/N_side_pts1;
h2 = abs(d-c)/N_side_pts2;

% Rectangular Grid Generation
side_pts1 = linspace(a,b,N_side_pts1)'* ones(1,N_side_pts2);
side_pts2 = [linspace(c,d,N_side_pts2)'*ones(1,N_side_pts1)]';

points= [];
for i = 1:N_side_pts1
    for j = 1:N_side_pts2
        points = [points; [side_pts1(i,j), side_pts2(i,j)]];
    end
end

% Perturbation
% We define how many standard deviations from the mean is the box size in

% we only want to perturb points not on the boundary.
% find those points
Bnd = and(and(points(:,1) ~= a, points(:,1) ~= b), ...
    and(points(:,2) ~= c, points(:,2) ~= d));

% the appropraite dimension
nStdDev = 10;
xPert = normrnd(0, h1/nStdDev, length(points(:,1)), 1).*Bnd;
yPert = normrnd(0, h2/nStdDev, length(points(:,1)), 1).*Bnd;

points = [points(:,1) + xPert, points(:,2) + yPert];

% remove points that are outside the boundary
in_pts = ~or(or(points(:,1) < a, points(:,1) > b),...
    or(points(:,2) < c, points(:,2) > d ));
points = points(find(in_pts),:);


% Visualize plot
scatter(points(:,1), points(:,2), 'ko')


%{
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

% Calculate unit normal vectors for each voronoi line

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
%}

