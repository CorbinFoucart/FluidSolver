% Using tic, toc to find matrix inversion times in Matlab

clear all; close all; clc; 

t = zeros(1,1000);
for n = 1:1000
    A = rand(n,n);
    b = rand(n,1);
    tic;
    x = A\b;
    t(n) = toc;
end
plot(t)