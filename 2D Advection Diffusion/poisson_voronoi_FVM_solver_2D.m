% Poisson FVM Solver, 2D
% Corbin Foucart
% ------------------------------------------------------------------------
% In this program, we read in the data from the previously created voronoi
% mesh and sovle the Poisson Equation in 2 dimensions over the mesh.
% Voronoi nodes are labeled 1:N. We look to solve AX = B.

clear all; close all; clc; 
format long 

% Read in Data, initialize column vectors and matrices. 
vfilename =  'voronoi_sqrgrid_data_30pts.mat';
load(vfilename);

normsave_filename = 'norm_data_sqrgrid';

N = length(points(:,1)); % size of matrix and column vectors
A = zeros(N);
B = zeros(N,1);

% Diffusive Coefficient
D = 1;

% ------------------- Intialize Boundary Conditions -------------------- %
% This is arbitrary. Chosen to check laplace equation first.
% Top and bottom of square 0, sides 10
% Amounts to finding sides of the square and recording those values as 10
% Also compute source integration terms at barycenter
for i = 1:N
    if dtri_points(i,3) == 1 % if boundary point
        %if dtri_points(i,1) == a ||  dtri_points(i,1) == b
        %    B(i) = 1;
        %end
        A(i,i) = 1;
    else % not boundary point
        %source_z = poisson_source_2d(points(i,1), points(i,2), D);
        source_z = poisson_source_2d(vertex_barycenters(i,1), vertex_barycenters(i,2), D);
        boxA = voronoi_areas(i);
        box_quad = boxA * source_z;
        B(i) = box_quad;
    end
end

% ------------------ Discretize Inner Voronoi Points ------------------- %
for i = 1:N
    if dtri_points(i,3) == 0 % if inner point
        nb_ind = find(neighbor_mat(i,:)); % indices of neighboring vertices
        for j = 1:length(nb_ind)
            vert_pts = [dtri_points(i,1:2); dtri_points(nb_ind(j),1:2)];
            dist = sqrt((vert_pts(1,1)-vert_pts(2,1))^2 + (vert_pts(1,2)-vert_pts(2,2))^2);
            % find edge length
            templns = edge_matrix(edge_matrix(:,1)== min(i, nb_ind(j)), :);
            for k = 1:length(templns(:,1))
                if templns(k,2) == max(i, nb_ind(j));
                    len = templns(k,3);               
                end
            end
            if len > 1e-12
                A(i,i) = A(i,i) + D*len/dist;
                A(i,nb_ind(j)) = A(i,nb_ind(j)) - D*len/dist;
            end
        end
    end
end

X = A\B;

solns = [dtri_points(:,1:2), X];
x = solns(:,1);
y = solns(:,2);
z = solns(:,3);

for i = 1:length(z)
    if abs(z(i)) < 1e-12
        z(i) = 0;
    end
end

% ------------------------- Error Measurement -------------------------- %

% L^2 norm
c_check = source_check(x,y);
diff = abs(c_check-z);
l2_norm = sqrt(sum(voronoi_areas.*(diff).^2));
m = matfile(normsave_filename, 'Writable', true);
m.l2_data = [m.l2_data; [length(points(:,1)), l2_norm]];
m.l2_data = unique(sortrows(m.l2_data, 1), 'rows');

% L infinity norm
linf_norm = max(diff);
m.linf_data = [m.linf_data; [length(points(:,1)), linf_norm]];
m.linf_data = unique(sortrows(m.linf_data, 1), 'rows');

% H1 norm
h1_norm = sqrt(diff'*A*diff);
m.h1_data = [m.h1_data; [length(points(:,1)), h1_norm]];
m.h1_data = unique(sortrows(m.h1_data, 1), 'rows');

mdat = load(normsave_filename)

sum_areas = sum(voronoi_areas)

% ---------------------- Plotting, Visualization  ---------------------- %

figure()
triplot = delaunay(x,y);
trisurf(triplot,x,y,z)
colormap hsv
alpha(.6)
