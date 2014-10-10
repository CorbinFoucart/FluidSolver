% Line segment, normal dot product quadrature
% Corbin Foucart
% -------------------------------------------------------------------

% takes as input:
% unit vector normal to segment
% vector values along the divided segment
% length of the line segment

function integral = adif2d_dotquad(uv, vvs, len)

% Computation
dotp = vvs(:,1)*uv(1) + vvs(:,2)*uv(2);

% Integration
% Riemmann
%integral = sum(dotp*ds)

% Trapezoidal Method
% Scaling by the length of the segment is handled in the here, not in the
% main function. 


integral = len*trapz(dotp);

% plotting, visualization, animation
%{
steps = length(vvs(:,1));
anim_time = 2; % seconds
time_step = anim_time/steps
max_vv = max(max(abs(vvs))); % For axis 

figure()
for i = 1:steps
    plot([0, vvs(i,1)],[0, vvs(i, 2)])
    hold all
    plot([0, uv(1)], [0, uv(2)], 'r-')
    plot([0, uv(1)*dotp(i)], [0, uv(2)*dotp(i)], 'ko')
    axis(1.1*[-max_vv, max_vv, -max_vv, max_vv])
    pause(time_step)
    clf()
end
close all;
%}

end