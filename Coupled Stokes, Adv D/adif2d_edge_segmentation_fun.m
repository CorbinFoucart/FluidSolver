% Advection Diffusion 2D solver: Edge Segmentation Function
% Corbin Foucart
% ------------------------------------------------------------------------
% Function that ...

% Output
% Cell array conatining 1x3 cells
% Each 'row' is each voronoi line
    % Within a 'row':
        % 'column' 1 is a matrix containing the points on the voronoi line
        % where intersections occur
        % 'column' 2 is a column array of the midpoints of the segments of
        % the voronoi line determined by the intersections
        % 'column' 3 is the triangle number that each segment determined by
        % its midpoint belongs to.

function out = adif2d_edge_segmentation_fun(v_lines, neighbor_mat,...
          edge_matrix, points, conn_list, adj_triangles) 


%clear all; close all; clc;
%Read in Data, initialize column vectors and matrices. 
%vfilename =  'voronoi_data_test.mat';
%load(vfilename);
        
for n = 1:length(v_lines(:,1))

% Build Collection of local points forming a polygon that encloses the
% endpoints of each voronoi line
vline_pts = [v_lines(n,1:2); v_lines(n,3:4)];
edge_pt_nums = edge_matrix(n,1:2);
edge_pts = points(edge_pt_nums, :);

% patch finder causing problems
[patch_vertices,  patch_tri_nums, patch_edge_nums] = ...
    adif2d_vline_patch_finder(vline_pts, edge_pt_nums, neighbor_mat,...
        conn_list, points, edge_matrix);
    
% Build point vectors to check for intersections
vline_int_pts = v_lines(n,1:4);
% edge points: select from edge matrix
edge_pt_nums = edge_matrix(patch_edge_nums,1:2);
edge_int_pts = zeros(length(edge_pt_nums(:,1)),4);
for i = 1:length(edge_pt_nums(:,1))
    edge_int_pts(i,:) = [points(edge_pt_nums(i,1),:),points(edge_pt_nums(i,2),:)];
end
   
out = lineSegmentIntersect(vline_int_pts, edge_int_pts);

% get indices in the adjacency matrix of intersecting lines
adj_ind = find(out.intAdjacencyMatrix);
% get corresponding edge numbers for those intersections
int_edge_nums = patch_edge_nums(adj_ind);

% case where v line is entirely within a triangle, find out which triangle
% the voronoi line is in, then we are done; must check all triangles
% because we don't have the edge numbers to guide us here. 
if isempty(adj_ind)
    midpt = 1/2*[sum(vline_pts(:,1)), sum(vline_pts(:,2))];
    dst = sqrt((vline_pts(1,1) - vline_pts(2,1))^2 ...
         + (vline_pts(2,2) - vline_pts(1,2))^2);
    for j = 1:length(patch_tri_nums)
        tri_pts = conn_list(patch_tri_nums(j),:);
            check_pts = points(tri_pts(1,:),:);
            [lm, tru] = adif2d_barycentric_coords(check_pts,...
                midpt);
            if tru == 1
                loc_triangle_num = patch_tri_nums(j);
                break
            end
    end
    crossing_data{n,1} = {[], midpt, dst, loc_triangle_num};
    
% Normal case, where we have some set of intersecting edges with the vline.
% Use the edges to get all adjacent triangles over the v line, and figure
% out which segment belongs to which. 
else
    seg_points = [vline_pts(1,:); [[out.intMatrixX(adj_ind)]',[out.intMatrixY(adj_ind)]'];...
        vline_pts(2,:)];
    % calculate distance from 1st endpoint
    seg_points(:,3) = sqrt((seg_points(:,1) - seg_points(1,1)).^2 +...
    (seg_points(:,2) - seg_points(1,2)).^2);
    % Append edge numbers for the intersections
    seg_points(:,4) = [0; int_edge_nums; 0];
    % Append the triangle numbers contained between the edges as the next
    % two columns
    seg_points(:,5:6) = [0,0; adj_triangles(int_edge_nums,:);0,0];
    % Sort the points by distance
    seg_points = sortrows(seg_points,3);

    % function that takes seg_points and returns the separate segments
    % along with the triangles that contain them

    % find midpoints of each segment
    seg_mid = [];
    dist = [];
    for k = 1:length(seg_points(:,1)) - 1
        seg_mid(k,:) = 1/2*[sum(seg_points(k:k+1, 1)), sum(seg_points(k:k+1, 2))];
        dist(k,1) = sqrt((seg_points(k+1,1) - seg_points(k,1))^2 ...
            + (seg_points(k+1,2) - seg_points(k,2))^2);
    end 

    % check 
    patch_tris = unique(adj_triangles(int_edge_nums,:));
    % remove zero entries from triangles on the boundary - there is no
    % triangle 0
    patch_tris(patch_tris == 0) = [];
    loc_triangle_nums = zeros(length(seg_mid(:,1)), 1);
    for i = 1:length(seg_mid(:,1))
        for j = 1:length(patch_tris)
            tri_pts = conn_list(patch_tris(j),:);
            check_pts = points(tri_pts(1,:),:);
            seg_mid(i,:);
            [lm, tru] = adif2d_barycentric_coords(check_pts,...
                seg_mid(i,:));
            if tru == 1
                loc_triangle_nums(i) = patch_tris(j);
                break
            end
        end
    end

    % Save data in cell
    % Each row is each voronoi line
    % Within a row:
        % column 1 is a matrix containing the points on the voronoi line
        % where intersections occur
        % column 2 is a column array of the midpoints of the segments of
        % the voronoi line determined by the intersections
        % column 3 is the triangle number that each segment determined by
        % its midpoint belongs to.
    crossing_data{n,1} = {seg_points, seg_mid, dist, loc_triangle_nums};
end

% for loop through voronoi edges end
end

out = crossing_data;
% function end
end