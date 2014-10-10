% Voronoi 2D mesh generator, adjacent triangle list
% Corbin Foucart
%-------------------------------------------------------------------------
% Input: edge matrix, connectivity list
% Output: rows: edge number, 
%         columns: nums of triangles adjacent to the edge

function return_adj_triangles = vg2d_adj_triangles_list_creation(...
    edge_matrix, conn_list)

adj_triangles = zeros(length(edge_matrix),2);
for i = 1:length(edge_matrix)
    % Get relevant vertices
    vertex1 = edge_matrix(i,1);
    vertex2 = edge_matrix(i,2);
    
    % Get boolean matrix of the connectivity list containing each vertex
    v1_search = (conn_list == vertex1);
    v2_search = (conn_list == vertex2);
    
    %  boolean matrix containing either value
    v_search = v1_search + v2_search;
    
    % sum along rows
    v_sum = sum(v_search, 2);
    
    % check for the rows containing 2s in the sum matrix
    % should ONLY be 1 or 2, for boundary edges or edges btwn triangles
    ind = [find(v_sum == 2)]';
    
    if length(ind) == 1 % boundary triangle
        adj_triangles(i,:) = [ind(1), 0];
    else adj_triangles(i,:) = ind;
    end 
end

return_adj_triangles = adj_triangles;

end

% Graveyard  Code, much slower method %
%{
adj_triangles = zeros(length(edge_matrix),2);
for i = 1:length(edge_matrix)
    candidate_triangles = [];
    epoints = edge_matrix(i,:);
    for k = 1:length(conn_list)
        if ismember(epoints(1), conn_list(k,:)) ...
                && ismember(epoints(2), conn_list(k,:))
            candidate_triangles = [candidate_triangles(1:end), k];
        end
    end
    if length(candidate_triangles) == 1
            candidate_triangles = [candidate_triangles(1:end), 0];
    end
    adj_triangles(i,:) = candidate_triangles ;
end
%}