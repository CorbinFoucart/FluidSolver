% L2 Norm Calculation
% Corbin Foucart
% -------------------------------------------------------------------------
% Takes in 2 vectors, the solution vector x in Ax=b, and the exact solution
% vector u.

function out = norm_l2_calc(x,u)

l2_norm = sqrt(sum((x-u).^2));
out = l2_norm;
end