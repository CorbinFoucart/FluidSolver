% Temperature interpolation Function
% Corbin Focuart

function T_l = temp_interp(Tsol, corner_nodes, corner_x, corner_y, xl, yl)

% lin interpolation of temperature using barycentric coords
R_bary = [corner_x, corner_y];
r_bary = [xl; yl];

% get temperature values at the corners of the triangle
corner_T = Tsol(corner_nodes);

lamb = adif2d_barycentric_coords(R_bary, r_bary);
T_l = lamb(1)*corner_T(1) + lamb(2)*corner_T(2) + lamb(3)*corner_T(3);  

end