% Save L2 norm function
% Corbin Foucart
% -------------------------------------------------------------------------
% Appends the result of a function to a previously stored matrix
% l2_data, linf_data, h1_data

function out = norm_save(N_points, norm, variableName)

% Read in Data
filename = 'norm_data'
m = matfile(filename, 'Writable', true)
m.variableName = [m.variableName; [N_points, norm]]
out = 1
end