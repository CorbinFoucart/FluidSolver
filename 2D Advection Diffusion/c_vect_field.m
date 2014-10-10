% Continuous Vector Field Function
% Corbin Foucart
% -------------------------------------------------------------------------
% Takes in (x,y) and returns vector r = [x,y]

function [rx ry] = c_vect_field(x,y)

%rx = ones(length(x),1);
%ry = 0.*y;

%rx = -x;
%ry = y;

% <u,v>
rx =  2.0.*(1.0 - x).^2.0.*x.^2.0.*(1.0 - y).^2.0.*y ...
     - 2.0.*(1.0 - x).^2.0.*x.^2.0.*(1.0 - y).*y.^2.0;

ry = -2.0.*(1.0 - x).^2.0.*x.*(1.0 - y).^2.0.*y.^2.0 ...
    + 2.*(1.0 - x).*x.^2.0.*(1.0 - y).^2.0.*y.^2.0;
end

