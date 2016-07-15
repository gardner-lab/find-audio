% this is a demo showing the use of our dynamic time warping package 
% we provide both Matlab version and C/MEX version
% the C/MEX version is much faster and highly recommended

% clear;clc;close all;

%compile;
mex CFLAGS='$FLAGS -Weverything' LDFLAGS='$LDFLAGS -framework Accelerate' dtw_path_c.c;

a = rand(10, 100);
b = rand(10, 120);

tic;
p1=dtw_path(a, b, 1.5);
t1=toc;

tic;
p2=dtw_path_c(a, b, 1.5);
t2=toc;

fprintf('Using Matlab version: distance=%d, running time=%f\n',size(p1, 2),t1);
fprintf('Using C/MEX version: distance=%d, running time=%f\n',size(p2, 2),t2);

% checks
disp(any(p1(:) ~= p2(:)));

% benchmark
tms = zeros(1, 1000);
for i = 1:1000
    tic;
    d2=dtw_ua_c(a,b);
    t2=toc;
    tms(i) = t2;
end
disp(mean(tms));
