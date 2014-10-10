% Norm plotter
% Corbin Foucart
% -------------------------------------------------------------------------
clear alll; close all; clc;

% Clearing mode

clearMode = 1;

% Initial norm data save;

if clearMode
    l2_data = []
    linf_data= []
    h1_data = []
    res_data= []
    filename = 'norm_data_sqrgrid'
    save(filename, 'h1_data', 'linf_data', 'l2_data', 'res_data')
else
    % Norm Plotting
    load('norm_data_sqrgrid', 'l2_data','linf_data', 'h1_data', 'res_data')
    h = 1./sqrt(l2_data(:,1));
    %figure()
    %loglog(h, l2_data(:,2))
    %title('mesh size h vs l2, log')
    %xlabel('mesh size h')
    %ylabel('L^2 norm')
    pl2 = polyfit(log(h), log(l2_data(:,2)),1)

    %figure()
    %loglog(h,linf_data(:,2))
    %title('mesh size h against L^\infty')
    %xlabel('mesh size h')
    %ylabel('L^\infty norm')
    plinf = polyfit(log(h), log(linf_data(:,2)), 1)

    %figure()
    %loglog(h,abs(h1_data(:,2)))
    %title('mesh size h against h1 norm')
    %xlabel('mesh size h')
    %ylabel('h1 norm')
    ph1 = polyfit(log(h), log(h1_data(:,2)), 1)
    
    res_data
    mean(res_data)
end

