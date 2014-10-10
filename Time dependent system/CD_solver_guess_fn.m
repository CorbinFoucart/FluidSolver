% Convection Diffusion FVM solver function
% Guess Mode
% Corbin Foucart
% ------------------------------------------------------------------------
% In this program, we read in the data from the previously created voronoi
% mesh and sovle the Convection Diffusion Eqation in 2 dimensions over the mesh.
% Voronoi nodes are labeled 1:N. We look to solve AX = B.

function [rT_vals, rA, rb] = ...
    CD_solver_guess_fn(vfilename)

% Read in mesh Data, initialize column vectors and matrices. 
load(vfilename);

D = 1; % diffusive constant

N = length(points(:,1)); % size of matrix and column vectors
A = zeros(N);
B = zeros(N,1);

% --------------------------- Guess Mode -------------------------------- %

% This is the case where we want to randomly generate temperature values as
% the first iteration. We also want to assemble a diffusion matrix, but not
% an convection matrix. 

% -------------------------- Vector Field Stuff ------------------------ %
% Get mesh crossing data
cross_dat = adif2d_edge_segmentation_fun(v_lines, neighbor_mat,...
            edge_matrix, points,conn_list, adj_triangles); 
        
% we know that this is we do not yet have velocity data
% from the crouzeix solver, and that we must use v = 0 everywhere for our
% initial guess.

% ------------------- Intialize Boundary Conditions -------------------- %
% Initialize BC, calcualte right hand side. 
for i = 1:N
    if dtri_points(i,3) == 1 % if boundary point
        if dtri_points(i,2) == c
            B(i) = 1;
        end
        A(i,i) = 1;
    else % not boundary point
        %source_z = zSource(vertex_barycenters(i,1), vertex_barycenters(i,2),D);
        %boxA = voronoi_areas(i);
        %box_quad = boxA * source_z;
        %B(i) = box_quad;
        
        B(i) = 0;
    end
end


% Diffusion
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
            % tolerance for when matlab assigns a finite value to the nodes
            % of zero length
            if len > 1e-12
                A(i,i) = A(i,i) + D*len/dist;
                A(i,nb_ind(j)) = A(i,nb_ind(j)) - D*len/dist;
            end
        end
    end
end

% Returning Temperature Guesses
guessVals = zeros(N,1);
for i = 1:N
    if dtri_points(i,3) == 1 % if boundary point
        if dtri_points(i,2) == c
            guessVals(i) = 1;
        end
    else % not boundary point
        guessVals(i) = rand(1,1);
    end
end

% function return values
rT_vals = guessVals;
rA = A;
rb = B;

end



