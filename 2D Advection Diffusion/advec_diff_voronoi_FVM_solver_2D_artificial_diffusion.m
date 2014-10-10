% Poisson FVM Solver, 2D
% Corbin Foucart
% ------------------------------------------------------------------------
% In this program, we read in the data from the previously created voronoi
% mesh and sovle the Poisson Equation in 2 dimensions over the mesh.
% Voronoi nodes are labeled 1:N. We look to solve AX = B.

clear all; close all; clc; 

% Read in Data, initialize column vectors and matrices. 
vfilename =  'voronoi_data_60pt';    %'voronoi_sqrgrid_data_5pts.mat' 'voronoi_data_5pt' 
load(vfilename);
normsave_filename = 'norm_data_sqrgrid';

D = 1; % diffusive constant

N = length(points(:,1)); % size of matrix and column vectors
A = zeros(N);
B = zeros(N,1);

% -------------------------- Vector Field Stuff ------------------------ %
% Get mesh crossing data
cross_dat = adif2d_edge_segmentation_fun(v_lines, neighbor_mat,...
            edge_matrix, points,conn_list, adj_triangles);
        
[vx, vy] = c_vect_field(points(:,1), points(:,2));

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
        source_z = poisson_source_2d(vertex_barycenters(i,1), vertex_barycenters(i,2),D);
        boxA = voronoi_areas(i);
        box_quad = boxA * source_z;
        B(i) = box_quad;
    end
end

% Advection
% for each vertex
z_param = zeros(length(v_lines(:,1)),1);
for j = 1:length(points(:,1))
if dtri_points(j,3) == 0 % if inner point
    
    % get voronoi lines
    loc_vlines = find(vertex_vlines(j,:));
    
    % for each voronoi line contained around each vertex
    for i = 1:length(loc_vlines)

        % get the vertex numbers to which the voronoi line correpsonds to
        c_vertex_numbers = edge_matrix(loc_vlines(i),1:2);
        
        % need to know which point is not the one we are currently dealing
        % with as a vertex in order to properly direction the unit vector.
        other_pt = c_vertex_numbers(find(c_vertex_numbers ~= j));
        
        % Get the actual points and distances
        c_vertex_points = points(c_vertex_numbers, :);
        dist = sqrt((c_vertex_points(2,1) - c_vertex_points(1,1))^2 ...
            + (c_vertex_points(2,2) - c_vertex_points(1,2))^2);

        % calculate the unit normal to the boundary
        uv = (points(other_pt,:) - points(j,:))/norm(points(other_pt,:) - points(j,:));
        
        c_seg = cross_dat{loc_vlines(i),1}{1,1};
        c_mids = cross_dat{loc_vlines(i),1}{1,2};
        c_dist = cross_dat{loc_vlines(i),1}{1,3};
        c_trinum = cross_dat{loc_vlines(i),1}{1,4};
        
        % for each segment of each voronoi line
        for k = 1:length(c_mids(:,1))
        % make sure voronoi segment has finite length
        if c_dist(k) ~= 0 && v_lines(loc_vlines(i),5) ~= 0
                % get segment endpoint vector values by barycentric
                % interpolation
                if isempty(c_seg)
                    % this is the case where the v line is totally included
                    % in one triangle, which is included in c_trinum
                    % endpoints of v_line should be interpolated over the
                    % triangle
                    lamb1 = adif2d_barycentric_coords(...
                    points(conn_list(c_trinum(k),:),:),...
                    v_lines(loc_vlines(i),1:2));
                    lamb2 = adif2d_barycentric_coords(...
                    points(conn_list(c_trinum(k),:),:), v_lines(loc_vlines(i),3:4));
                else
                lamb1 = adif2d_barycentric_coords(...
                    points(conn_list(c_trinum(k),:),:), c_seg(k,1:2));
                lamb2 = adif2d_barycentric_coords(...
                    points(conn_list(c_trinum(k),:),:), c_seg(k+1,1:2));
                end
                
                % get velocity vectors at the vertices of the triangle
               corn_v = [vx(conn_list(c_trinum(k),:)), vy(conn_list(c_trinum(k),:))];
               int_v1 = [dot(corn_v(:,1), lamb1), dot(corn_v(:,2), lamb1)];
               int_v2 = [dot(corn_v(:,1), lamb2), dot(corn_v(:,2), lamb2)];
               quadr = adif2d_dotquad(uv, [int_v1; int_v2], c_dist(k));
                    
               zp = abs(quadr./v_lines(loc_vlines(i),5).*edge_matrix(loc_vlines(i),3)./D);
               z_param(loc_vlines(i)) = z_param(loc_vlines(i)) + zp;
               
               % j is our vertex number, other_pt is the neighboring vertex
               % number
               A(j,j) = A(j,j) + 1/2*quadr;            
               A(j,other_pt) = A(j,other_pt) + 1/2*quadr;
               
        end
        end

    end
end
end

% Diffusion
% ------------------ Discretize Inner Voronoi Points ------------------- %
save_Diff = zeros(size(A));
eig_check = zeros(size(A));
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
            
            % find edge number 
            en = find((sum((edge_matrix(:,1:2) == i) + (edge_matrix(:,1:2) == nb_ind(j)), 2)) == 2);
            locz = z_param(en);
            % ask Dr. Linke
            if locz ~= 0
                art_diff_param = locz./2.*coth(locz./2) - 1;
            else
                art_diff_param = 0;
            end
            
            A(i,i) = A(i,i) + D*len/dist + D*art_diff_param*len/dist;
            A(i,nb_ind(j)) = A(i,nb_ind(j)) - D*len/dist - D*art_diff_param*len/dist;
            
            % For H1 norm, we must save ONLY the diffusion elements (not
            % the artificial diffusion!)
            save_Diff(i,i) = save_Diff(i,i) + D*len/dist;
            save_Diff(i,nb_ind(j)) = save_Diff(i, nb_ind(j)) - D*len/dist;
            
            % for checking eigenvalues, we need only the artificial
            % diffusion matrix
            eig_check(i,i) = eig_check(i,i) + D*art_diff_param*len/dist;
            eig_check(i,nb_ind(j)) = eig_check(i,nb_ind(j)) - D*art_diff_param*len/dist;
        end
    else
        eig_check(i,i) = D;
    end
end

X = A\B;

solns = [dtri_points(:,1:2), X];
x = solns(:,1);
y = solns(:,2);
z = solns(:,3);

diff = abs(z - source_check(x,y));

% ------------------------- Error Measurement -------------------------- %
%
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
h1_norm = sqrt(diff'*save_Diff*diff);
m.h1_data = [m.h1_data; [length(points(:,1)), h1_norm]];
m.h1_data = unique(sortrows(m.h1_data, 1), 'rows');

res = A * z - B;
maxerr = max(abs(res));
residual = 1./length(res)*sum(abs(res));

mdat = load(normsave_filename)

sum_areas = sum(voronoi_areas)
%
% ---------------------- Plotting, Visualization  ---------------------- %

%
figure()
triplot = delaunay(x,y);
trisurf(triplot,x,y,z)
colormap hsv
axis square
alpha(.6)
%




