% Barycenter Calculator Function
% Takes in set of vertices, orders them into a convex hull, and passes the
% points to another polygon geometry program. 

function centroid = barycenter_calc(x,y)

vi = convhull(x,y);
centroid = polygeom(x(vi), y(vi));

end