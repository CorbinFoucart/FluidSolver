% Poisson FVM Solver, 2D
% Corbin Foucart
% ------------------------------------------------------------------------
% In this program, we read in the data from the previously created voronoi
% mesh and sovle the Poisson Equation in 2 dimensions over the mesh.
% Voronoi nodes are labeled 1:N. We look to solve AX = B.

clear all; close all; clc; 

% Read in Data, initialize column vectors and matrices. 
vfilename =  'voronoi_data_20pt_un.mat';
load(vfilename);

N = length(points(:,1)); % size of matrix and column vectors
A = zeros(N);
B = zeros(N,1);

% ------------------- Intialize Boundary Conditions -------------------- %
% This is arbitrary. Chosen to check laplace equation first.
% Top and bottom of square 0, sides 10
% Amounts to finding sides of the square and recording those values as 10

for i = 1:N
    if dtri_points(i,3) == 1 % if boundary point
        if dtri_points(i,1) == a || dtri_points(i,1) == b
            B(i) = 1;
        end
        A(i,i) = 1;
    end
end

% ------------------ Discretize Inner Voronoi Points ------------------ %


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
            A(i,i) = A(i,i) - len/dist;
            A(i,nb_ind(j)) = A(i,nb_ind(j)) + len/dist;
        end
    end
end

X = A\B;

solns = [dtri_points(:,1:2), X]
x = solns(:,1);
y = solns(:,2);
z = solns(:,3);

% ---------------------- Plotting, Visualization  ---------------------- %


figure()
triplot = delaunay(x,y);
trisurf(triplot,x,y,z)
colormap hsv
alpha(.6)

