% crouzeix.m in 2D

clear all; close all; clc;
% output format
% format longEng
% format long

% options
isRT0mode = 1;

% for visualization
h1 = figure;

% machine eps
cmpEps = 1.0e-12;

% definition of the academic velocity solution
u = @(x, y) 2 * (1 - x)^2 * x^2 * (1 - y)^2 * y - 2 * (1 - x)^2 * x^2 * (1 - y) * y^2;
v = @(x, y) -2 * (1 - x)^2 * x * (1 - y)^2 * y^2 + 2 * (1 - x) * x^2 * (1 - y)^2 * y^2;
ux = @(x,y) 4 * x * (1 - 3 * x + 2 * x^2) * y * (1 - 3 * y + 2 * y^2);
uy = @(x, y) 2 * (-1 + x)^2 * x^2 * (1 - 6 * y + 6 * y^2);
vx = @(x, y) -2 * (1 - 6 * x + 6 * x^2) * (-1 + y)^2 * y^2;
vy = @(x, y) -4 * x * (1 - 3 * x + 2 * x^2) * y * (1 - 3 * y + 2 * y^2);

% definition of the academic pressure solution
p = @(x, y) 0;

% definition of the right hand side
f1 = @(x, y) 12 * x.^2 - 24 * x.^3 + 12 * x.^4 - 4 * y + 24 * x * y - 48 *x.^2 * y + 48 * x.^3 * y ... 
 - 24 * x.^4 * y + 12 * y.^2 - 72 * x * y.^2 + 72 * x.^2 * y.^2 - 8 * y.^3 + 48 * x * y.^3 ... 
 - 48 * x.^2 * y.^3;
f2 = @(x, y) 4 * x - 12 * x.^2 + 8 * x.^3 - 24 * x * y + 72 * x.^2 * y - 48 * x.^3 * y - 12 * y.^2 + 48 * x * y.^2 - 72 * x.^2 * y.^2 + 48 * x.^3 * y.^2 + 24 * y.^3 - 48 * x * y.^3 - 12 * y.^4 + 24 * x * y.^4;

% definition of the nodes in the mesh
%nx = [0 1 0 1]'; % [0 1 0 1]'; % [0 1 0 2 2]'; % [0 1 0 2]';
%ny = [0 0 1 1]'; % [0 0 1 1]'; % [0 0 1 2 1]'; % [0 0 1 2]';
%nx = [0 1/2 1  0  1/2  1  0 1/2 1]';
%ny = [0  0  0 1/2 1/2 1/2 1  1  1]';

%nx = [0  0  1/2 3/10 35/100 7/10 65/100 1  1  0 1/2 1]';
%ny = [0 1/2  0  3/10 65/100 7/10 35/100 0 1/2 1  1  1]';

[nx, ny] = getStructuredSquareGrid(9);

nodenum = length(nx);

% mesh generation
trirep = DelaunayTri(nx, ny);

e = trirep.edges();

edgenum = length(e);

trinum = trirep.size(1);

simplex2edges = zeros(trinum, 3);
hi = zeros(trinum, 1);

rhs = zeros(2 * edgenum + trinum, 1);

% defining the boundary
isBoundaryEdge = zeros(edgenum, 1);
for i=1:edgenum
    xa = nx(e(i, 1));
    ya = ny(e(i, 1));
    
    xb = nx(e(i, 2));
    yb = ny(e(i, 2));
    
    if ((abs(xa) < cmpEps) && (abs(xb) < cmpEps)) ...
            || ((abs(ya) < cmpEps) && (abs(yb) < cmpEps)) ...
            || ((abs(xa - 1.0) < cmpEps) && (abs(xb - 1.0) < cmpEps)) ...
            | ((abs(ya - 1.0) < cmpEps) && (abs(yb - 1.0) < cmpEps))
        isBoundaryEdge(i) = 1;
    end
end

% definition of quadrature formulas

% order 1
%qmat = [1/3 1/3 1/3];
%ws = [1]';

% order 2
%qmat = [1/2 1/2 0; 0 1/2 1/2; 1/2 0 1/2];
%ws = [1/3 1/3 1/3]';

% order 3
%qmat = [1/3 1/3 1/3; 0.6 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.6];
%ws = [-27/48 25/48 25/48 25/48];

% order 10
%qmat = [1/3 1/3 1/3; ...
%    0.4978654329544748 0.4978654329544748 0.00426913409105033; ...
%    0.4978654329544748 0.00426913409105033 0.4978654329544748; ...
%    0.00426913409105033 0.4978654329544748 0.4978654329544748; ...
%    0.4280124497290562 0.4280124497290562 0.14397510054188760; ...
%    0.4280124497290562 0.14397510054188760 0.4280124497290562; ...
%    0.14397510054188760 0.4280124497290562 0.4280124497290562; ...
%   0.1847564127432246 0.1847564127432246 0.6304871745135509; ...
%    0.1847564127432246 0.6304871745135509 0.1847564127432246; ...
%    0.6304871745135509 0.1847564127432246 0.1847564127432246; ...
%   0.0204812185716776 0.0204812185716776 0.9590375628566449; ...
%    0.0204812185716776 0.9590375628566449 0.0204812185716776; ...
%    0.9590375628566449 0.0204812185716776 0.0204812185716776; ...
%    0.8284234338466945 0.0350029898972720 0.1365735762560335; ...
%    0.8284234338466945 0.1365735762560335 0.0350029898972720; ...
%    0.1365735762560335 0.8284234338466945 0.0350029898972720; ...
%    0.1365735762560335 0.0350029898972720 0.8284234338466945; ...
%    0.0350029898972720 0.8284234338466945 0.1365735762560335; ...
%    0.0350029898972720 0.1365735762560335 0.8284234338466945; ...
%    0.6297073291529187 0.3327436005886387 0.037549070258443; ...
%    0.6297073291529187 0.037549070258443 0.3327436005886387; ...
%    0.3327436005886387 0.6297073291529187 0.037549070258443; ...
%    0.3327436005886387 0.037549070258443 0.6297073291529187; ...
%    0.037549070258443 0.6297073291529187 0.3327436005886387; ...
%    0.037549070258443 0.3327436005886387 0.6297073291529187; ...
%    ];
%ws = [0.08352339980519638 ...
%    0.007229850592056743 0.007229850592056743 0.007229850592056743 ...
%    0.07449217792098051 0.07449217792098051 0.07449217792098051 ...
%    0.07864647340310853 0.07864647340310853 0.07864647340310853 ...
%    0.006928323087107504 0.006928323087107504 0.006928323087107504 ...
%    0.02951832033477940 0.02951832033477940 0.02951832033477940 ...
%    0.02951832033477940 0.02951832033477940 0.02951832033477940 ...
%    0.03957936719606124 0.03957936719606124 0.03957936719606124 ...
%    0.03957936719606124 0.03957936719606124 0.03957936719606124 ...
%    ];

% order 15
qmat = [ 0.08343840726174993 0.4582807963691250 0.4582807963691250; ...
    0.4582807963691250 0.4582807963691250 0.08343840726174993; ...
    0.4582807963691250 0.08343840726174993 0.4582807963691250; ...
    0.19277907084173887 0.4036104645791306 0.4036104645791306; ...
    0.4036104645791306 0.4036104645791306 0.19277907084173887; ...
    0.4036104645791306 0.19277907084173887 0.4036104645791306; ...
    0.2931971679130254 0.2931971679130254 0.4136056641739493; ...
    0.2931971679130254 0.4136056641739493 0.2931971679130254; ...
    0.4136056641739493 0.2931971679130254 0.2931971679130254; ...
    0.1464677869427729 0.1464677869427729 0.7070644261144541; ...
    0.1464677869427729 0.7070644261144541 0.1464677869427729; ...
    0.7070644261144541 0.1464677869427729 0.1464677869427729; ...
    0.0563628676656034 0.0563628676656034 0.8872742646687931; ...
    0.0563628676656034 0.8872742646687931 0.0563628676656034; ...
    0.8872742646687931 0.0563628676656034 0.0563628676656034; ...
    0.0165751268583703 0.0165751268583703 0.9668497462832593; ...
    0.0165751268583703 0.9668497462832593 0.0165751268583703; ...
    0.9668497462832593 0.0165751268583703 0.0165751268583703; ...
    0.009912203309225 0.2395345541547944 0.7505532425359808; ...
    0.009912203309225 0.7505532425359808 0.2395345541547944; ...
    0.2395345541547944 0.009912203309225 0.7505532425359808; ...
    0.2395345541547944 0.7505532425359808 0.009912203309225; ...
    0.7505532425359808 0.009912203309225 0.2395345541547944; ...
    0.7505532425359808 0.2395345541547944 0.009912203309225; ...
    0.015803770630228 0.4048788073183400 0.5793174220514320; ...
    0.015803770630228 0.5793174220514320 0.4048788073183400; ...
    0.4048788073183400 0.015803770630228 0.5793174220514320; ...
    0.4048788073183400 0.5793174220514320 0.015803770630228; ...
    0.5793174220514320 0.015803770630228 0.4048788073183400; ...
    0.5793174220514320 0.4048788073183400 0.015803770630228; ...
    0.005143608816971 0.0950021131130449 0.8998542780699844; ...
    0.005143608816971 0.8998542780699844 0.0950021131130449; ...
    0.0950021131130449 0.005143608816971 0.8998542780699844; ...
    0.0950021131130449 0.8998542780699844 0.005143608816971; ...
    0.8998542780699844 0.005143608816971 0.0950021131130449; ...
    0.8998542780699844 0.0950021131130449 0.005143608816971; ...
    0.0489223257529888 0.1497531073222740 0.8013245669247372; ...
    0.0489223257529888 0.8013245669247372 0.1497531073222740; ...
    0.1497531073222740 0.0489223257529888 0.8013245669247372; ...
    0.1497531073222740 0.8013245669247372 0.0489223257529888; ...
    0.8013245669247372 0.0489223257529888 0.1497531073222740; ...
    0.8013245669247372 0.1497531073222740 0.0489223257529888; ...
    0.0687687486325192 0.2869196124413350 0.6443116389261458; ...
    0.0687687486325192 0.6443116389261458 0.2869196124413350; ...
    0.2869196124413350 0.0687687486325192 0.6443116389261458; ...
    0.2869196124413350 0.6443116389261458 0.0687687486325192; ...
    0.6443116389261458 0.0687687486325192 0.2869196124413350; ...
    0.6443116389261458 0.2869196124413350 0.0687687486325192; ...
    0.1684044181246992 0.2818356680990846 0.5497599137762162; ...
    0.1684044181246992 0.5497599137762162 0.2818356680990846; ...
    0.2818356680990846 0.1684044181246992 0.5497599137762162; ...
    0.2818356680990846 0.5497599137762162 0.1684044181246992; ...
    0.5497599137762162 0.1684044181246992 0.2818356680990846; ...
    0.5497599137762162 0.2818356680990846 0.1684044181246992; ...
    ];
ws = [ 0.03266181884880529 0.03266181884880529 0.03266181884880529...
    0.02741281803136436 0.02741281803136436 0.02741281803136436 ...
    0.02651003659870330 0.02651003659870330 0.02651003659870330 ...
    0.02921596213648611 0.02921596213648611 0.02921596213648611 ...
    0.01058460806624399 0.01058460806624399 0.01058460806624399 ...
    0.003614643064092035 0.003614643064092035 0.003614643064092035 ...
    0.008527748101709436 0.008527748101709436 0.008527748101709436 ...
    0.008527748101709436 0.008527748101709436 0.008527748101709436 ...
    0.01391617651669193 0.01391617651669193 0.01391617651669193 ...
    0.01391617651669193 0.01391617651669193 0.01391617651669193 ...
    0.004291932940734835 0.004291932940734835 0.004291932940734835 ...
    0.004291932940734835 0.004291932940734835 0.004291932940734835 ...
    0.01623532928177489 0.01623532928177489 0.01623532928177489 ...
    0.01623532928177489 0.01623532928177489 0.01623532928177489 ...
    0.02560734092126239 0.02560734092126239 0.02560734092126239 ...
    0.02560734092126239 0.02560734092126239 0.02560734092126239 ...
    0.03308819553164567 0.03308819553164567 0.03308819553164567 ...
    0.03308819553164567 0.03308819553164567 0.03308819553164567 ...
    ];

pdelibqmat = zeros(54, 2);
for i=1:54
    pdelibqmat(i, 1) = -qmat(i, 1) + qmat(i, 2) - qmat(i, 3);
    pdelibqmat(i, 2) = -qmat(i, 1) - qmat(i, 2) + qmat(i, 3);
end

% construct adjacency matrix returning for each triangle
% all adjacent edges
for i=1:edgenum
    tris = edgeAttachments(trirep, e(i, 1), e(i, 2));
    trisvec=cell2mat(tris);
    for j=1:length(trisvec)
        k=1;
        t = trisvec(j);
        while(k > 0 && k <= hi(t))
            if simplex2edges(t, k+1) == i
              k = -1;
            else
              k = k+1;
            end
        end
        if k > 0
            hi(t) = hi(t) + 1;
            simplex2edges(t, hi(t)) = i;
        end
    end
end

simplex2nodes = zeros(trinum, 3);
hi = zeros(trinum, 1);

% construct adjacency matrix returning for each triangle
% all adjacent nodes
for i=1:nodenum
    tris = vertexAttachments(trirep, i);
    trisvec=cell2mat(tris);
    for j=1:length(trisvec)
        k=1;
        t = trisvec(j);
        while(k > 0 && k <= hi(t))
            if simplex2nodes(t, k+1) == i
              k = -1;
            else
              k = k+1;
            end
        end
        if k > 0
            hi(t) = hi(t) + 1;
            simplex2nodes(t, hi(t)) = i;
        end
    end
end

% compute triangle areas by Heron's formula
triareas = zeros(trinum, 1);
for i=1:trinum
    n1 = simplex2nodes(i, 1);
    n2 = simplex2nodes(i, 2);
    n3 = simplex2nodes(i, 3);
    
    l12 = sqrt((nx(n1) - nx(n2)).^2 + (ny(n1) - ny(n2)).^2);
    l13 = sqrt((nx(n1) - nx(n3)).^2 + (ny(n1) - ny(n3)).^2);
    l23 = sqrt((nx(n2) - nx(n3)).^2 + (ny(n2) - ny(n3)).^2);
    
    s = (l12 + l13 + l23) / 2;
    
    triareas(i) = sqrt(s * (s - l12) * (s - l13) * (s - l23));
end

xbcs = zeros(trinum, 1);
ybcs = zeros(trinum, 1);
% for RT0 reconstruction and pressure visualization
for i=1:trinum
    n1 = simplex2nodes(i, 1);
    n2 = simplex2nodes(i, 2);
    n3 = simplex2nodes(i, 3);
    
    xbcs(i) = 1/3 * (nx(n1) + nx(n2) + nx(n3));
    ybcs(i) = 1/3 * (ny(n1) + ny(n2) + ny(n3));
end

%simplex2edges
%simplex2nodes
%triareas

triplot(simplex2nodes, nx, ny);

% allocating the stiffness matrix
is = [];
js = [];
ss = [];
mat = sparse(is, js, ss, 2 * edgenum + trinum, 2 * edgenum + trinum);

phi_as = zeros(3, 1);
phi_bs = zeros(3, 1);
phi_cs = zeros(3, 1);

phi_mat = zeros(3);
phi_rhs = zeros(3, 1);
phi_sol = zeros(3, 1);

for i=1:trinum
    e1 = simplex2edges(i, 1);
    e2 = simplex2edges(i, 2);
    e3 = simplex2edges(i, 3);
    
    emp1x = (nx(e(e1, 1)) + nx(e(e1, 2))) / 2;
    emp1y = (ny(e(e1, 1)) + ny(e(e1, 2))) / 2;
    
    emp2x = (nx(e(e2, 1)) + nx(e(e2, 2))) / 2;
    emp2y = (ny(e(e2, 1)) + ny(e(e2, 2))) / 2;
    
    emp3x = (nx(e(e3, 1)) + nx(e(e3, 2))) / 2;
    emp3y = (ny(e(e3, 1)) + ny(e(e3, 2))) / 2;
    
    phi_mat(1, 1) = 1;
    phi_mat(2, 1) = 1;
    phi_mat(3, 1) = 1;
    phi_mat(1, 2) = emp1x;
    phi_mat(2, 2) = emp2x;
    phi_mat(3, 2) = emp3x;
    
    phi_mat(1, 3) = emp1y;
    phi_mat(2, 3) = emp2y;
    phi_mat(3, 3) = emp3y;
    
    phi_sol = phi_mat \ [1 0 0]';
    phi_as(1) = phi_sol(1);
    phi_bs(1) = phi_sol(2);
    phi_cs(1) = phi_sol(3);
    
    phi_sol = phi_mat \ [0 1 0]';
    phi_as(2) = phi_sol(1);
    phi_bs(2) = phi_sol(2);
    phi_cs(2) = phi_sol(3);
    
    phi_sol = phi_mat \ [0 0 1]';
    phi_as(3) = phi_sol(1);
    phi_bs(3) = phi_sol(2);
    phi_cs(3) = phi_sol(3);
    
    %phi_as
    %phi_bs
    %phi_cs
    
    le = [e1 e2 e3]';
    
    % assembly of the Laplacian operators
    for ii=1:3
        for jj=1:3
            if(~isBoundaryEdge(le(ii)) && ~isBoundaryEdge(le(jj)))
              mat(le(jj), le(ii)) = mat(le(jj), le(ii)) + triareas(i) * (phi_bs(ii) * phi_bs(jj) + phi_cs(ii) * phi_cs(jj));
              mat(edgenum+le(jj), edgenum+le(ii)) = mat(edgenum+le(jj), edgenum+le(ii)) + triareas(i) * (phi_bs(ii) * phi_bs(jj) + phi_cs(ii) * phi_cs(jj));
            end
        end
    end
    
    % assembly of the divergence and the pressure gradient operators
    for ii=1:3
        % (-u_x, q)
        if (~isBoundaryEdge(le(ii))) && (i ~= 1)
          mat(2 * edgenum + i, le(ii)) = mat(2 * edgenum + i, le(ii)) + triareas(i) * (-phi_bs(ii));
          mat(le(ii), 2 * edgenum + i) = mat(le(ii), 2 * edgenum + i) + triareas(i) * (-phi_bs(ii));
        end
        
        % (-v_y, q)
        if (~isBoundaryEdge(le(ii))) && (i ~= 1)
          mat(2 * edgenum + i, edgenum + le(ii)) = mat(2 * edgenum + i, edgenum + le(ii)) + triareas(i) * (-phi_cs(ii));
          mat(edgenum + le(ii), 2 * edgenum + i) = mat(edgenum + le(ii), 2 * edgenum + i) + triareas(i) * (-phi_cs(ii));
        end
    end
    
    % assembly of the right hand side
    [qN, dummy] = size(qmat);
    
    if ~isRT0mode  
      for qi=1:qN
          xc = baryToCart(trirep, [i]', [qmat(qi, 1) qmat(qi, 2) qmat(qi, 3)]);
          xl = xc(1);
          yl = xc(2);
        
          for ii=1:3
            if ~isBoundaryEdge(le(ii))
              rhs(le(ii)) = rhs(le(ii)) + triareas(i) * ws(qi) * f1(xl, yl) * (phi_as(ii) + xl * phi_bs(ii) + yl * phi_cs(ii));
              rhs(edgenum + le(ii)) = rhs(edgenum + le(ii)) + triareas(i) * ws(qi) * f2(xl, yl) * (phi_as(ii) + xl * phi_bs(ii) + yl * phi_cs(ii));
            end
          end
      end
    end % if ~isRT0mode
    
    if isRT0mode
        % RT0 reconstruction of the form v = a_T + c_T ( x - x_T)
        for qi=1:qN
            xc = baryToCart(trirep, [i]', [qmat(qi, 1) qmat(qi, 2) qmat(qi, 3)]);
            xl = xc(1);
            yl = xc(2);
            
            for ii=1:3
                if ~isBoundaryEdge(le(ii))
                    xa = nx(e(le(ii), 1));
                    ya = ny(e(le(ii), 1));
                    
                    xb = nx(e(le(ii), 2));
                    yb = ny(e(le(ii), 2));
                    
                    % get edge midpoint
                    empx = 1/2 * (xa + xb);
                    empy = 1/2 * (ya + yb);
                    
                    % compute triangle edge length sigma
                    sigma = sqrt((xb-xa)^2 + (yb-ya)^2);
                    
                    % b parameter Raviart
                    b = sigma/triareas(i);
                    
                    % compute unit normal to triangle edge sigma:
                    nsigmax = -(yb - ya) / sigma;
                    nsigmay = (xb - xa) / sigma;
                    
                    % get right orientation of edge normal
                    if (empx - xbcs(i)) * nsigmax ...
                            + (empy - ybcs(i)) * nsigmay < 0
                        nsigmax = -nsigmax;
                        nsigmay = -nsigmay;
                    end
                    
                    % nsigmax
                    % nsigmay
                    
                    % horizontal degree of freedom
                    aTx = 1/triareas(i) * sigma * nsigmax * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * nsigmax * (empy - ybcs(i));
                    
                    % cT is 1/2 of the divergence of hat function (0, 1)
                    cT = phi_bs(ii) / 2;
                    
                    rhs(le(ii)) = rhs(le(ii)) + triareas(i) * ws(qi) ...
                        * (f1(xl, yl) * (aTx + cT * (xl - xbcs(i))) ...
                    + f2(xl, yl) * (aTy + cT * (yl - ybcs(i))));
                
                    % vertical degree of freedom
                    aTx = 1/triareas(i) * sigma * nsigmay * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * nsigmay * (empy - ybcs(i));
                    
                    % cT is 1/2 of the divergence of hat function (0, 1)
                    cT = phi_cs(ii) / 2;
                    rhs(edgenum + le(ii)) = rhs(edgenum + le(ii)) + triareas(i) * ws(qi) ...
                        * (f1(xl, yl) * (aTx + cT * (xl - xbcs(i))) ...
                    + f2(xl, yl) * (aTy + cT * (yl - ybcs(i))));
                
                    % Construct Raviart Thomas Array for each edge
                    % Row contains [aTx aTy b unormx unormy]
                    
                    aTx = 1/triareas(i) * sigma * (empx - xbcs(i));
                    aTy = 1/triareas(i) * sigma * (empy - ybcs(i));
                    
                    RT{i,1}{ii} = [aTx, aTy, b, nsigmax, nsigmay, le(ii)];
                else
                    RT{i,1}{ii} = [0, 0, 0, 0, 0, le(ii)];
                end
            end
        end
    end % if isRT0mode
end

% assembly of the boundary behavior
for i=1:edgenum
    if(isBoundaryEdge(i))
        mat(i, i) = 1;
        mat(edgenum + i, edgenum + i) = 1;
    end
end
mat(2 * edgenum + 1, 2 * edgenum + 1) = -1;

% linear solve
sol = mat \ rhs;

% flow visualization
xmps = zeros(edgenum, 1);
ymps = zeros(edgenum, 1);
usol = zeros(edgenum, 1);
vsol = zeros(edgenum, 1);
psol = zeros(trinum,  1);

for i=1:edgenum
    xmps(i) = 1/2 * (nx(e(i, 1)) + nx(e(i,2)));
    ymps(i) = 1/2 * (ny(e(i, 1)) + ny(e(i,2)));
    usol(i) = sol(i);
    vsol(i) = sol(edgenum + i);
end

% raviart thomas velocity computation
RT_Ax = zeros(length(triareas),1);
RT_Ay = zeros(length(triareas),1);
RT_b  = zeros(length(triareas),1);
for i = 1:length(triareas)
    for j = 1:3
        localv = [usol(RT{i,1}{1,j}(6)), vsol(RT{i,1}{1,j}(6))];
        % second row of cell relevant velocities
        RT{i,1}{2,j} = localv;
        % 3rd row dot product between normals and velocities
        RT{i,1}{3,j} = dot(RT{i,1}{1,j}(4:5), localv);
        % 4th row dot products multiplied with the respective factors
        % scalar * array
        RT{i,1}{4,j} = RT{i,1}{3,j}*RT{i,1}{1,j}(1:3);
    end
    % next row in cell the values [aTx aTy b] for each triangle
    RT{i,2} = RT{i,1}{4,1} + RT{i,1}{4,2} + RT{i,1}{4,3}; 
    
    % Save values to more easily accessible vectors
    RT_Ax(i,1) =  RT{i,2}(1);
    RT_Ay(i,1) =  RT{i,2}(2);
    RT_b(i,1)  =  RT{i,2}(3);
end


% change pressure level
pavg = 0;
mesomega = 0;
for i=1:trinum
    mesomega = mesomega + triareas(i);
    pavg = pavg + triareas(i) * sol(2 * edgenum + i);
end
pavg = pavg / mesomega;

for i=1:trinum
    psol(i) = sol(2 * edgenum + i) - pavg;
end



% compute norms
l2normExact = 0;
l2norm = 0;
l2errnorm = 0;

h1normExact = 0;
h1norm = 0;
h1errnorm = 0;

divnorm = 0;

pnormExact = 0;
pnorm = 0;
perrnorm = 0;

for i=1:trinum
    e1 = simplex2edges(i, 1);
    e2 = simplex2edges(i, 2);
    e3 = simplex2edges(i, 3);
    
    emp1x = (nx(e(e1, 1)) + nx(e(e1, 2))) / 2;
    emp1y = (ny(e(e1, 1)) + ny(e(e1, 2))) / 2;
    
    emp2x = (nx(e(e2, 1)) + nx(e(e2, 2))) / 2;
    emp2y = (ny(e(e2, 1)) + ny(e(e2, 2))) / 2;
    
    emp3x = (nx(e(e3, 1)) + nx(e(e3, 2))) / 2;
    emp3y = (ny(e(e3, 1)) + ny(e(e3, 2))) / 2;
    
    phi_mat(1, 1) = 1;
    phi_mat(2, 1) = 1;
    phi_mat(3, 1) = 1;
    phi_mat(1, 2) = emp1x;
    phi_mat(2, 2) = emp2x;
    phi_mat(3, 2) = emp3x;
    
    phi_mat(1, 3) = emp1y;
    phi_mat(2, 3) = emp2y;
    phi_mat(3, 3) = emp3y;
    
    phi_sol = phi_mat \ [1 0 0]';
    phi_as(1) = phi_sol(1);
    phi_bs(1) = phi_sol(2);
    phi_cs(1) = phi_sol(3);
    
    phi_sol = phi_mat \ [0 1 0]';
    phi_as(2) = phi_sol(1);
    phi_bs(2) = phi_sol(2);
    phi_cs(2) = phi_sol(3);
    
    phi_sol = phi_mat \ [0 0 1]';
    phi_as(3) = phi_sol(1);
    phi_bs(3) = phi_sol(2);
    phi_cs(3) = phi_sol(3);
    
    le = [e1 e2 e3]';
    
    % evaluate function u(x, y) and v(x, y) at quadrature points
    for qi=1:qN
        xc = baryToCart(trirep, [i]', [qmat(qi, 1) qmat(qi, 2) qmat(qi, 3)]);
        xl = xc(1);
        yl = xc(2);
        
        l2normExact = l2normExact + triareas(i) * ws(qi) * (u(xl, yl)^2 + v(xl, yl)^2);
        h1normExact = h1normExact + triareas(i) * ws(qi) ...
            * (ux(xl, yl)^2 + uy(xl, yl)^2 + vx(xl, yl)^2 + vy(xl, yl)^2);
        pnormExact = pnormExact + triareas(i) * ws(qi) * p(xl, yl)^2;
        
        uu = 0;
        vv = 0;
        
        gux = 0;
        guy = 0;
        gvx = 0;
        gvy = 0;
        for ii=1:3
            uu = uu + usol(le(ii)) * (phi_as(ii) + xl * phi_bs(ii) + yl * phi_cs(ii));
            vv = vv + vsol(le(ii)) * (phi_as(ii) + xl * phi_bs(ii) + yl * phi_cs(ii));
            
            gux = gux + usol(le(ii)) * phi_bs(ii);
            guy = guy + usol(le(ii)) * phi_cs(ii);
            gvx = gvx + vsol(le(ii)) * phi_bs(ii);
            gvy = gvy + vsol(le(ii)) * phi_cs(ii);
        end
        l2norm = l2norm + triareas(i) * ws(qi) * (uu^2 + vv^2);
        l2errnorm = l2errnorm + triareas(i) * ws(qi) * ((u(xl, yl) - uu)^2 + (v(xl, yl) - vv)^2);
        
        h1norm = h1norm + triareas(i) * ws(qi) * (gux^2 + guy^2 + gvx^2 + gvy^2);
        h1errnorm = h1errnorm + triareas(i) * ws(qi) ...
            * ((ux(xl, yl) - gux)^2 + (uy(xl, yl) - guy)^2 ...
             + (vx(xl, yl) - gvx)^2 + (vy(xl, yl) - gvy)^2);
         
        divnorm = divnorm + triareas(i) * ws(qi) * (gux + gvy)^2;
        pnorm = pnorm + triareas(i) * ws(qi) * (psol(i))^2;
        perrnorm = perrnorm + triareas(i) * ws(qi) * (p(xl, yl) - psol(i))^2;
    end
end

save('RAx10', 'RT_Ax')
save('RAy10', 'RT_Ay')

%{
l2normExact = sqrt(l2normExact)
l2norm = sqrt(l2norm)
l2errnorm = sqrt(l2errnorm)

h1normExact = sqrt(h1normExact)
h1norm = sqrt(h1norm)
h1errnorm = sqrt(h1errnorm)

divnorm = sqrt(divnorm)

pnormExact = sqrt(pnormExact)
pnorm = sqrt(pnorm)
perrnorm = sqrt(perrnorm)
%}



h2 = figure;
quiver(xmps, ymps, usol, vsol);

figure()
quiver(xbcs, ybcs, RT_Ax, RT_Ay)

h3 = figure;

ptri = delaunay(xbcs, ybcs);
trisurf(ptri, xbcs, ybcs, psol);







