% Voronoi Patch Finder
% Corbin Foucart
% -----------------------------------------------------------------------
% Finds a small region of points/triangles/edges around a specified
% voronoi line

function [rn_patch_vertices, rn_tri_nums, rn_edge_nums] = ...
    adif2d_vline_patch_finder(vline_pts, edge_pt_nums, neighbor_mat,...
        conn_list, points, edge_matrix)
    

while(true)
    % take vertex nums in array edge_pt_nums and get all the
    % neighboring points to the vertex set
    neigh_nums = [];
    for i = 1:length(edge_pt_nums)
        neigh_nums = unique([neigh_nums, find(neighbor_mat(edge_pt_nums(i),:))]);
    end
    % point numbers in the region
    neigh_nums = neigh_nums';
    
    % Unnecessary?
    % actual points in the region
    %local_pts = points(neigh_nums,:);

    % Pass points to a function that returns all the d. triangle
    % numbers conatined in the region. 
    local_tris = adif2d_tri_returner(neigh_nums, conn_list);

    % check each triangle to see if each voronoi line endpoint is
    % contained in one of the returned triangles

    truth = zeros(2,1);
    for i = 1:2
        for j = 1:length(local_tris)
            % extract point numbers conn_list
            tst_pt_nums = conn_list(local_tris, :);
            % get actual points for each row
            test_acpts = points(tst_pt_nums(j,:),:);
            % calcualte barycentric coords, and thus whether point is
            % inside
            [lamb, in] = adif2d_barycentric_coords(test_acpts, vline_pts(i,:));
            % if endpoint located in one of the triangles, note and
            % start looking at the other endpoint
            if in == 1
                truth(i) = 1;
                break
            end
        end
    end

    % If endpoints contained within our triangle patch, we're done. If they
    % don't, we need to take the neighbor list as the new edge_pts to get
    % the neighbors of the neighbors on the next iteration. 
    if sum(truth) == 2
        break
    else edge_pt_nums = neigh_nums;
    end
end

rn_edge_nums = adif2d_edge_returner(neigh_nums, edge_matrix);
rn_patch_vertices = neigh_nums;
rn_tri_nums = local_tris;

end