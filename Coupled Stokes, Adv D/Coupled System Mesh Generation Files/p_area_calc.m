
function p_area = p_area_calc(x,y)

vi = convhull(x,y);
p_area = polyarea(x(vi), y(vi));

end