% online plotter

clear all;
close all;
clc;

cd 'RT1 - data files'

filename = 'Ra50000mesh20temp.mat'
tempSave = strcat(filename(1:end - 4),'.jpg')
velSave = strcat(filename(1:end - 8),'vel.jpg')
load(filename)

%cd 'RT1 - Images'


% triangulations for plotting
triplot = delaunay(x,y);

% Gridfitting for countour plots
contNodesX = 0:0.1:1;
contNodesY = 0:0.1:1;

% magnitude of vectors
%contZ = sqrtm(usol*usol' + vsol*vsol');
%{
contZ = usol + vsol;

[contx, conty, contz] = gridfit(xmps, ymps, contZ,...
    contNodesX, contNodesY);
%}
% Temperature Profile
figure()
trisurf(triplot,x,y,Tsol)
title('Temperature Profile')
colormap jet
axis square
alpha(.7)
xlabel('x')
ylabel('y')
zlabel('T')

h = gcf;
saveas(h,tempSave);

% Velocity Field
figure()
quiver(xmps, ymps, usol, vsol);
title('Velocity Profile')
axis square
alpha(.6)
xlabel('x')
ylabel('y')

h = gcf
saveas(h,velSave)

%{
% Contour Plot of Velocity Field
figure()
contour(contx, conty, contz);
title('Velocity Contour Profile')
axis square
xlabel('x')
ylabel('y')
%}

cd ..

